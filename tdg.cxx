#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <deque>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cerrno>
#include <cassert>
#include <unistd.h>
#include "matrix.hxx"
#include "mpi_buffer.hxx"
#include "streams.hxx"


#define CHECKMPI(fn, str) \
  { int const ret = (fn); if (ret != MPI_SUCCESS) mpierror (ret, (str)); }


#define TYPE_MSG 'M' //! A message. See MSG_*. 
#define TYPE_STKELEM 'S' //! Stack element, as response to a request. 
#define TYPE_TOKEN  'T' //! Token for ADUV.
#define TYPE_DONOR 'D' //! Answer to a donor request.

#define MSG_EOC 'E' //! End of computation. 
#define MSG_REQ 'R' //! Request work. 
#define MSG_DENY 'D' //! Deny work. 
#define MSG_DREQ 'O' //! Request donor from P1. 

#define TOKEN_BLACK 'B'
#define TOKEN_WHITE 'W'
#define TOKEN_NONE 'N'

//#define TAG_CAN_WAIT 1
#define TAG_NEEDS_ATTENTION 20


/* Message receive buffer. */
size_t const BUFSIZE = 0xffffu;
void * recv_buf = 0;
size_t recv_buf_len;
/* MPI status. */
MPI_Status status;

//! GCZ-AHD counter. 
unsigned donor;
//! ADUV 
char mycolor, token;
//! "Would give out work" flag. 
bool wouldgive = true;
//! Time */
double t1 = 0, t2 = 0;
//! MPI rank of process.
int rank = -1;
//! Size of the world.
int worldsize = -1;
//! Minimal unit of work done before probing for messages.
size_t const BATCHSIZE = 1000u;
//! Batch counter.
unsigned batch_size = BATCHSIZE;

#ifdef NDEBUG
null_ostream log;
#else
prefixed_ostream log (std::cerr);
#endif // NDEBUG


/*! 
  \brief Structure representing one element of algorithm's DFS stack.  */
struct StackElem
{
  //! Representation of X and Y sets of edges.
  std::vector<bool> colour;
  //! Bitmap whose bits indicate validity of bits in above bitmap.
  std::vector<bool> valid;
  /*! 
    Offset of the rightmost 1 in the future new element generated from this
    element.  */
  unsigned next;
  //! Number of edges from all valid edges that have been used so far.
  unsigned used_edges;

  StackElem ()
    : next (std::numeric_limits<unsigned>::max ()), 
      used_edges (std::numeric_limits<unsigned>::max ())
  { }

  StackElem (unsigned width, unsigned nxt, unsigned used)
    : colour (width), valid (width), next (nxt), used_edges (used)
  { }

  bool operator == (StackElem const & x) const
  {
    return colour == x.colour && valid == x.valid;
  }

  bool operator < (StackElem const & x) const
  {
    if (valid < x.valid)
      return true;
    else if (valid == x.valid)
      if (colour < x.colour)
        return true;
      else
        return false;
    else
      return false;
  }
};


/**
   \brief Serializes stack element.
 */
void
serialize_stackelem (mpi_obuffer & obuf, StackElem const & el)
{
  /*
    Layout: |width|next|used_edges|edge_1|colour_1|...
    ...|edge_{used_edges}|colour_{used_edges}|
   */
  size_t const sz = el.colour.size ();
  obuf << sz << el.next << el.used_edges;
  for (unsigned i = 0, cnt = 0; i < sz && cnt < el.used_edges; ++i)
    {
      if (el.valid[i])
        {
          obuf << i << static_cast<char>(el.colour[i]);
          cnt += 1;
        }
    }
}


void
deserialize_stackelem (StackElem & el, mpi_ibuffer & ibuf)
{
  size_t sz;
  unsigned next, used_edges;
  ibuf >> sz >> next >> used_edges;
  el.colour.resize (sz);
  el.valid.resize (sz);
  el.next = next;
  el.used_edges = used_edges;
  for (unsigned i = 0; i < used_edges; ++i)
    {
      unsigned edge;
      char colour;
      ibuf >> edge >> colour;
      el.colour[edge] = colour;
      el.valid[edge] = true;
    }
}


/*!
  \brief Structure that holds elements of edge -> (vertex, vertex) mapping.  */
struct EdgeInfo
{
  //! First vertex of an edge.
  unsigned v1;
  //! Second vertex of an edge.
  unsigned v2;
  //! Next valid vertex in edgeinfo vector.
  unsigned next;

  EdgeInfo ()
    : v1 (std::numeric_limits<unsigned>::max ()), 
      v2 (std::numeric_limits<unsigned>::max ()), 
      next (std::numeric_limits<unsigned>::max ())
  { }
  
  EdgeInfo (unsigned iv1, unsigned iv2, unsigned nxt = 0)
    : v1 (iv1), v2 (iv2), next (nxt)
  { }
};


//! Graph type.
typedef UTMatrix<bool> Graph;
//! Set of vertices.
typedef std::set<unsigned> SetOfVertices;
//! Vector of sets of vertices' neighbours.
typedef std::vector<SetOfVertices> Neighbours;
//! Vector of bits that indicate presence of some edge in the graph.
typedef std::vector<bool> ValidEdges;
//! Mapping of pairs of vertices to sets of their common successors.
typedef std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> > PairsToCommonSucc;


//! Number of vertices in graph.
unsigned N = 0;
//! Maximal number of edges in graph.
unsigned max_edges = 0;
//! Number of edges in graph.
unsigned true_edge_count = 0;
//! Number of edges that we will use for colouring.
unsigned edge_count = 0;
//! First edge.
unsigned first_edge = 0;
//! Number of generated states.
unsigned states = 1;
//! Upper triangle matrix representing the graph.
Graph * graph = 0;

/*! 
  Vector of sets of neighbouring vertices for each vertex. For fast triangle
  testing.  */
Neighbours neighbours;

//! Map of pairs of vertices to their common successors.
PairsToCommonSucc common_succ;

/*!
  \brief Bitmap mapping integers for all possible edges to booleans,
  edges[num] == true iff the edge is present in the graph.  */
ValidEdges edges;
//! Uninteresting edges that are not part of any triangle.
ValidEdges unin_edges;

/*!
  \brief Bitmap with bits set to 1 for edges that we don't consider during
  colouring because they can't influence the outcome, i.e. edges that do not
  take part in any triangle.  */
ValidEdges precoloured_edges;

/*!
  \brief Map from edge numbers to pairs of vertices. edgeinfo[num] is only
  valid iff edges[num] == true.  */
std::vector<EdgeInfo> edgeinfo;
//! DFS algoritm's stack.
std::deque<StackElem> stack;
//! Element of result.
StackElem * result_element = 0;


/**
   Prints msg and possible message for errno to stderr and
   exits with EXIT_FAILURE.
   @param msg user supplied message
*/
void
error (const char * msg)
{
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if (! errno)
    log << msg << std::endl;
  else
    log << msg << ": " << strerror (errno) << std::endl;
  MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
  abort (); /* Just to silenece "warning: `noreturn' function does return" */
}


/**
   Prints msg and possible message for MPI error to stderr and
   exits with EXIT_FAILURE.
   @param msg user supplied message
*/
void
mpierror (int ret, const char * msg)
{
  int len, rank;
  char str[1000];

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Error_string (ret, str, &len);
  log << msg << ": " << str << std::endl;
  MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
  abort ();  /* Just to silenece "warning: `noreturn' function does return" */
}


/*!
  \brief Return index of element for element at position (x, y) if it is
  stored in upper triangular matrix of size NxN.  */
inline unsigned
index_for_xy (unsigned x, unsigned y)
{
  if (x < 1 || x > N || y < 1 || y > N)
    abort ();

  if (y > x)
    std::swap (x, y);
  return ((y - 1) * (N + (N - (y-1-1)))) / 2
    + (x - (y - 1)) - 1;
}


/*!
   \brief Initializes variables.  */
void 
init_structures ()
{
  // Graph.
  graph = new Graph (N);
  // Maximal number of edges in the graph.
  max_edges = (N * (N + 1ul)) / 2;
  // Vector of bools indicating valid edges.
  edges.resize (max_edges);
  unin_edges.resize (max_edges);
  // Map from edges to vertices.
  edgeinfo.resize (max_edges);
  // Vertices info.
  neighbours.resize (N);
}


/*!
  \brief Checks whether adding edge results into creating a monochromatic
  triangle.

  \param el StackElem to which we want to add edge.
  \param edge The edge being added.
  \return Pair of bools. When first == true then adding the edge would result
  in a triangle of the first colour. Similarly with the second colour.  */
std::pair<bool, bool>
check_triangle (StackElem const & el, unsigned edge)
{
  assert (edges[edge]);

  unsigned const v1 = edgeinfo[edge].v1;
  unsigned const v2 = edgeinfo[edge].v2;
  // Loops cannot create triangles with other edges.
  if (v1 == v2)
    return std::make_pair (false, false);
  
  // Find set of common successors.
  PairsToCommonSucc::const_iterator set_it = 
    common_succ.find (std::make_pair (v1, v2));
  assert (set_it != common_succ.end ());
  std::vector<unsigned>::const_iterator it = set_it->second.begin ();
  std::vector<unsigned>::const_iterator const it_end = set_it->second.end ();
  
  /*
    Test if candidate edges have the same colour. Loop until we either test
    all vertices or until we determine that adding vertex of either colour
    would create a triangle.  */
  std::pair<bool, bool> res (false, false);
  for (;
       it != it_end && !(res.first && res.second);
       ++it)
    {
      unsigned const u = *it;
      unsigned const e1 = index_for_xy (v1, u);
      unsigned const e2 = index_for_xy (v2, u);
      
      // Are both candidate edges in current set of added and coloured edges?
      if (! el.valid[e1] || ! el.valid[e2])
        continue;

      // Are both edges the same colour?
      if (el.colour[e1] == el.colour[e2])
        // Yes, what colour would be the triangle?
        if (el.colour[e1] == false)
          res.first = true;
        else
          res.second = true;
    }

  return res;
}


/*!
   \brief Read graph from file and initializes all necessary variables.  
   
   \param infile Input stream.  */
void
read_graph (std::istream & infile)
{
  // Read graph's size from the input file.
  infile >> N;
  // Allocate and initialize structures.
  init_structures ();
  // Read the graph from input the file.
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= N; ++j)
      {
        // Read edge from file.
        unsigned val;
        infile >> val;

        // Ignore everything that is bellow matrix' diagonal. 
        if (i > j)
          continue;

        // Update graph.
        (*graph)(j, i) = val;
        // Is there an edge here?
        if (val)
          {
            // Update number of edges.
            true_edge_count += 1;

            // Record the neighbourness of the two vertices.
            neighbours[j-1].insert (i);
            neighbours[i-1].insert (j);
            
            // Update bitmap of valid edges.
            unsigned const idx = index_for_xy (j, i);
            edges[idx] = true;
          }
      }
  std::cout << std::endl;
  
  // Construct edgeinfo.
  unsigned idx;
  unsigned last = std::numeric_limits<unsigned>::max ();
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= N; ++j)
      {
        // Ignore everything that is bellow matrix' diagonal. 
        if (i > j || !(*graph)(j, i))
          continue;
        
        idx = index_for_xy (j, i);
        edgeinfo[idx].v1 = j;
        edgeinfo[idx].v2 = i;

        // Compute common successors for both vertices of an edge.
        std::vector<unsigned> common;
        std::set_intersection 
          (neighbours[i-1].begin (), neighbours[i-1].end (),
           neighbours[j-1].begin (), neighbours[j-1].end (),
           std::back_inserter (common));

        if (common.size () == 0 || j == i)
          {
            // Mark uninteresting edge and continue.
            unin_edges[idx] = true;
            edges [idx] = false;
            continue;
          }
        edge_count += 1;

        // Insert new set of vertices into cache.
        common_succ.insert (std::make_pair (std::make_pair (j, i), common));
        
        // Is this the first edge that we see?
        if (last == std::numeric_limits<unsigned>::max ())
          {
            // Yes. No successor edge yet, next loop will set it.
            edgeinfo[idx] = EdgeInfo (j, i);
            first_edge = idx;
          }
        else
          {
            // Set successor of the previous edge.
            edgeinfo[last].next = idx;
            // Successor of this edge will be set in the next loop.
            edgeinfo[idx] = EdgeInfo (j, i);
            
            // Debug.
            if (rank == 0)
              std::cout << " successor edge " << idx 
                        << (unin_edges[last] ? " not interesting" : "")
                        << std::endl;
          }
        
        // Debug.
        if (rank == 0)
          std::cout << "edge " << idx 
                    << " vertices (" << j << ", " << i << ")";

        // Record this edge for the next loop.
        last = idx;
      }
  // The last edge doesn't have any successor.
  if (last != std::numeric_limits<unsigned>::max ())
    edgeinfo[last].next = std::numeric_limits<unsigned>::max ();
  
  if (rank == 0)
    {
      // Print uninteresting vertices/edge.
      std::cout << std::endl << "Edges that are not part of any triangle: ";
      for (unsigned i = 1; i <= N; ++i)
        for (unsigned j = 1; j <= N; ++j)
          {
            // Ignore everything that is bellow matrix' diagonal. 
            if (i > j || !(*graph)(j, i))
              continue;
            unsigned const idx = index_for_xy (i, j);
            if (unin_edges[idx])
              std::cout << idx << "/(" << i << "," << j << ") ";
          }
      std::cout << std::endl;
      
      // Some statistics.
      std::cout << std::endl;
      std::cout << "N: " << N << std::endl;
      std::cout << "edges: " << true_edge_count << std::endl;
      std::cout << "interesting edges: " << edge_count << std::endl;
      std::cout << "uninteresting edges: " << true_edge_count - edge_count 
                << std::endl;
    }
}


/*!
   \brief Prints the solution represented by stack element.  

   \param solution Stack element representing the solution.  */
void
print_solution (StackElem & solution)
{
  //std::cout << std::endl << "Solution has been found:" << std::endl;
  std::cout << std::endl << "[" << rank << "] Solution has been found:" 
            << std::endl;

  // Merge solution with pre-coloured edges.
  for (size_t i = 0; i < max_edges; ++i)
    {
      if (unin_edges[i])
        solution.valid[i] = true;
    }

  // Output first set of edges.
  std::cout << "red edges:";
  for (unsigned i = 0; i < max_edges; ++i)
    if (solution.valid[i] && solution.colour[i] == false)
      std::cout << " " << i;
  std::cout << std::endl;

  // Output second set of edges.
  std::cout << "white edges:";
  for (unsigned i = 0; i < max_edges; ++i)
    if (solution.valid[i] && solution.colour[i] == true)
      std::cout << " " << i;
  std::cout << std::endl;

  // Some statistics.
  std::cout << std::endl;
  std::cout << "generated states: " << states << std::endl;
}


/**
   MPI initialization.
*/
void
initialize_mpi (int & argc, char **& argv, int & rank, int & size)
{
  CHECKMPI (MPI_Init (&argc, &argv), "MPI_Init");
  CHECKMPI (MPI_Comm_rank (MPI_COMM_WORLD, &rank), "MPI_Comm_rank");
#ifndef NDEBUG
  {
    std::ostringstream oss;
    oss << "[" << rank << "] ";
    log.set_prefix (oss.str ());
  }
#endif
  log << "reporting" << std::endl;
  CHECKMPI (MPI_Comm_size(MPI_COMM_WORLD, &size), "MPI_Comm_size");

  // Allocate some memory for the incoming messages.
  recv_buf_len = BUFSIZE;
  recv_buf = new char [recv_buf_len];
}


void
initialize ()
{
  //! Receive buffer for MPI?

  if (rank == 0)
    token = TOKEN_WHITE;
  else
    token = TOKEN_NONE;
  mycolor = TOKEN_WHITE;
  if (rank == 0)
    {
      log << "Size of the world: " << worldsize << std::endl;
      // Put initial state on the stack.
      stack.push_back (StackElem (max_edges, first_edge, 0));
    }
}


void
end_computation (bool got_result)
{
  if (got_result)
    {
      t2 = MPI_Wtime ();
        // Prepare the message.
      mpi_obuffer obuf;
      obuf << TYPE_MSG << MSG_EOC;
      log << "sending MSG_EOC:" << std::endl;
      // Send end of computation.
      for (int i = 0; i < worldsize; ++i)
        {
          if (i == rank)
            continue;

          CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, i, 
                              TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), 
                    "MPI_Send");
          log << "\tnode " << i << std::endl;
        }
      // Print out the solution. 
      print_solution (*result_element);
      // Print the time.
      std::cout.setf (std::ios::fixed);
      std::cout << std::endl << "TIME: " << (t2 - t1) << std::endl;
    }
  
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize ();
  exit (EXIT_SUCCESS);
}


void
do_tokens ()
{
  log << "taking care of tokens" << std::endl;
  mpi_obuffer obuf;
  if (rank == 0)
    {
      obuf << TYPE_TOKEN << TOKEN_WHITE;
      log << "sending WHITE token to " << ((rank + 1) % worldsize) 
          << std::endl;
      CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, 
                          (rank + 1) % worldsize, TAG_NEEDS_ATTENTION, 
                          MPI_COMM_WORLD), "MPI_Send");
      token = TOKEN_NONE;
    }
  else
    {
      if (token != TOKEN_NONE)
        {
          obuf << TYPE_TOKEN << token;
          log << "sending '" << token << "' token to " 
              << ((rank + 1) % worldsize) << std::endl;
          CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED,
                              (rank + 1) % worldsize, TAG_NEEDS_ATTENTION,
                              MPI_COMM_WORLD), "MPI_Send");
          mycolor = TOKEN_WHITE;
          token = TOKEN_NONE;
        }
    }
}


void
deny_work_request (int from)
{
  // Prepare the message.
  mpi_obuffer obuf;
  obuf << TYPE_MSG << MSG_DENY;
  // Send it.
  CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, from, 
                      TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), "MPI_Send");
}


/**
   \brief Tries to find suitable donor stack element.
*/
StackElem *
find_suitable_elem ()
{
  unsigned const cut_height = edge_count * 2 / 3;
  StackElem * el;
  unsigned el_index;
  
  // Find suitable element from that we can donate work.
  if (stack.size () > cut_height)
    {
      el = &stack[cut_height];
      el_index = cut_height;
    }
  else
    {
      el = &stack.back ();
      el_index = stack.size () - 1;
    }
  // Check that the element can be expanded.
  if (el->next < max_edges)
    return el;
  
  /*
    Nope. We have to try harder to find suitable element.  Try to find one in
    direction of the botton of the stack.  */
  std::deque<StackElem>::iterator it = stack.begin ();
  std::advance (it, el_index);
  for (std::deque<StackElem>::reverse_iterator rit (it); 
       rit != stack.rend (); ++rit)
    {
      if (rit->next < max_edges)
        return &*rit;
    }

  /*
    Previous loop hasn't found any suitable element. Try find one in the
    direction of top of the stack.  */
  it = stack.begin ();
  std::advance (it, el_index);
  for (; it != stack.end (); ++it)
    {
      if (it->next < max_edges)
        return &*it;
    }

  error ("find_suitable_elem: Shouldn't happen.");
  abort ();
}


/* Either send out work or deny the request. */
void
process_work_request (int from)
{
  /*
     Do we have anything to give?
     Do we want to give at all?
  */
  if (stack.empty () || ! wouldgive)
    {
      /* Nope, deny the request. */
      if (! wouldgive)
        log << "we don't want to give anything, denying request" << std::endl;
      else
        log << "there is nothing to give, denying request" << std::endl;
      deny_work_request (from);
      return;
    }

  /* We have something to give. */
  StackElem * const el = find_suitable_elem ();
  unsigned const half = 2 * (edge_count - el->used_edges) / 2;
  if (half != 0 && rank > from)
    /* Change token. */
    mycolor = TOKEN_BLACK;
  /* Generate the half. */
  log << "generating " << half << " new stack elements from element that has"
      << " " << el->used_edges << " used edges (out of " << edge_count << ")"
      << std::endl;
  unsigned cnt = 0;
  std::list<StackElem> tmplist;
  while (!(cnt > half || el->next >= max_edges))
    {
      // Is it possible to go deeper in DFS tree?
      if (el->next >= max_edges)
        // No, return instead.
        break;
      // Save the number of the soon to be used edge.
      unsigned const inserted = el->next;
      // Adjust the pointer to next edge.
      el->next = edgeinfo[inserted].next;

      /* el:    [1 0 0 ... 0]
         |
         v
         newel: [1 1 0 ... 0]  */
      StackElem newel (*el);
      // Increment used edges counter for the new element.
      newel.used_edges = el->used_edges + 1;
      
      /*
        Check whether the new state is a valid one, i.e. that adding the edge
        with either colour won't result in a monochromatic triangle.  */
      std::pair<bool, bool> const res = check_triangle (newel, inserted);

      /*
        Push newel with one colour onto DFS stack if it doesn't create a
        monochromatic triangle.  */
      if (res.first == false)
        {
          newel.colour[inserted] = false;
          newel.valid[inserted] = true;
          tmplist.push_back (newel);
          states += 1;
          batch_size += 1;
          cnt += 1;
        }

      /*
        Push newel with the other colour onto DFS stack if it doesn't create
        monochromatic triangle.  */
      if (res.second == false)
        {
          newel.colour[inserted] = true;
          newel.valid[inserted] = true;
          tmplist.push_back (newel);
          states += 1;
          batch_size += 1;
          cnt += 1;
        }
    }
  log << "generated " << cnt << " elements" << std::endl;
  // Deny the request if we were not able to generate valid states.
  if (cnt == 0)
    {
      deny_work_request (from);
      return;
    }
  
  /* Send the half to requester. */
  mpi_obuffer obuf;
  // Message type and number of elements.
  obuf << TYPE_STKELEM << cnt;
  for (std::list<StackElem>::reverse_iterator it = tmplist.rbegin ();
       it != tmplist.rend (); ++it)
    {
      serialize_stackelem (obuf, *it);
    }
  // Finally the send.
  CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, from, 
                      TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), "MPI_Send");
}


void
process_donor_request (unsigned from)
{
  if (rank == 0)
    {
      mpi_obuffer obuf;
      obuf << TYPE_DONOR << donor;
      log << "request for donor has been received"
          << ", sending donor " << donor << " to process " << from 
          << std::endl;
      donor = (donor + 1) % worldsize;
      CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, from, 
                          TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), "MPI_Send");
      return;
    }
  else
    error ("Donor request to process != P1!!!");
}



void process_priority_message (mpi_ibuffer &);

unsigned
request_donor (void)
{
  /* Send request for donor to 0. */
  mpi_obuffer obuf;
  obuf << TYPE_MSG << MSG_DREQ;
  log << "sending request for donor to 0" << std::endl;
  CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, 0, 
                      TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), "MPI_Send");
  /* Wait for answer from 0. */
  mpi_ibuffer ibuf;
  while (1)
    {
      CHECKMPI (MPI_Recv (recv_buf, recv_buf_len, MPI_PACKED, MPI_ANY_SOURCE,
                          TAG_NEEDS_ATTENTION, MPI_COMM_WORLD, &status),
                "MPI_Recv");
      ibuf.attach (recv_buf, recv_buf_len);
      char type;
      ibuf >> type;
      if (type != TYPE_DONOR || status.MPI_SOURCE != 0)
        /* Some other priority message. */
        {
          ibuf.reset ();
          process_priority_message (ibuf);
          continue;
        }
      else
        /* Received the donor. */
        break;
    }

  char dnr;
  ibuf >> dnr;
  return dnr;
}


void
receive_stackelem (mpi_ibuffer & ibuf)
{
  unsigned cnt;
  ibuf >> cnt;
  log << "\treceived " << cnt << " elements" << std::endl;
  for (unsigned i = 0; i < cnt; ++i)
    {
      StackElem newel;
      deserialize_stackelem (newel, ibuf);
      stack.push_back (newel);
    }
}


/* Processes messages with TAG_NEEDS_ATTENTION. */
void
process_priority_message (mpi_ibuffer & ibuf)
{
  char type;
  char msg_type;

  ibuf >> type;
  log << "processing message type='" << type << "'" << std::endl;
  switch (type)
    {
    case TYPE_MSG:
      ibuf >> msg_type;
      log << "processing simple message type='" << msg_type << "'" 
          << std::endl;
      switch (msg_type)
        {
        case MSG_REQ:
          process_work_request (status.MPI_SOURCE);
          return;

        case MSG_DREQ:
          process_donor_request (status.MPI_SOURCE);
          return;

        case MSG_EOC:
          log << "end of computation has been received" << std::endl;
          //! \todo Add MPI_Barier for sync?
          end_computation (false);
          //MPI_Finalize ();
          //exit (EXIT_SUCCESS);

        default:
          log << "unhandled priority message!!!" << std::endl;
          error ("Unhandled MSG_ in process_priority_message()!!!");
        }

    case TYPE_TOKEN:
      {
        char tok;

        ibuf >> tok;
        log << "received '" << tok << "' token" << std::endl;
        if (rank == 0)
          if (tok == TOKEN_WHITE)
            {
              log << "got WHITE token back, invoking end_computation()" 
                  << std::endl;
              end_computation (false);
            }
          else
            {
              log << "coloring to 'W'" << std::endl;
              token = TOKEN_WHITE;
            }
        else
          if (mycolor == TOKEN_WHITE)
            {
              log << "color='" << mycolor << "' '" << token 
                  << "' has been received" << std::endl;
              token = tok;
            }
          else
            {
              log << "color='" << mycolor << "', coloring token to '"
                  << TOKEN_BLACK << "'" << std::endl;
              token = TOKEN_BLACK;
            }
        return;
      }
      
    case TYPE_STKELEM:
      receive_stackelem (ibuf);
      return;

    default:
      error ("Unhandled TYPE_* in process_priority_message()!!!");
    }
}


void 
request_work (int from)
{
  /* Send the request. */
  mpi_obuffer obuf;
  obuf << TYPE_MSG << MSG_REQ;
  log << "sending request for work to " << from << std::endl;
  CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, from, 
                      TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), "MPI_Send");
  
  /* Process answer. */
  while (1)
    {
      char type, msg_type;
      char dnr;

      CHECKMPI (MPI_Recv (recv_buf, recv_buf_len, MPI_PACKED, MPI_ANY_SOURCE,
                          TAG_NEEDS_ATTENTION, MPI_COMM_WORLD, &status),
                "MPI_Recv");
      mpi_ibuffer ibuf (recv_buf, recv_buf_len);
      if (status.MPI_SOURCE != from)
        {
          process_priority_message (ibuf);
          continue;
        }
      ibuf >> type;
      switch (type)
        {
        case TYPE_MSG:
          ibuf >> msg_type;
          switch (msg_type)
            {
              /*
            case MSG_EOE:
              log << "received MSG_EOE from " << from << "n";
              return;
              */

            case MSG_DENY:
              log << "received denying answer from " << from << std::endl;
              do_tokens ();
              do
                {
                  if (rank != 0)
                    /* Request a donor from 0. */
                    dnr = request_donor ();
                  else
                    {
                      dnr = donor;
                      donor = (donor + 1) % worldsize;
                    }
                  if (dnr != rank)
                    /* Send request to obtained donor and read results. */
                    request_work (dnr);
                  //continue;
                  /* Request a new donor from 0. */
                  //dnr = request_donor ();
                  /* Send request to obtained donor and read results. */
                  //request_work (dnr);
                }
              while (dnr == rank);
              return;

            default:
              log << "processing priority message from " << status.MPI_SOURCE
                  << " in request_work()" << std::endl;
              ibuf.reset ();
              process_priority_message (ibuf);
              continue;
            }

        case TYPE_STKELEM:
          {
            log << "received TYPE_STKELEM" << std::endl;
            receive_stackelem (ibuf);
            return;
          }

        default:
          log << "processing priority message from " << status.MPI_SOURCE
              << " in request_work()" << std::endl;
          ibuf.reset ();
          process_priority_message (ibuf);
          continue;
        }
    }
  error ("BUG!!! You should not have ever seen this!!!");
}
            
            
int
main (int argc, char * argv[])
try
{
  initialize_mpi (argc, argv, rank, worldsize);
  
  // Some basic checks and initialization.
  if (argc != 2)
    error ("Syntax: tdg <input graph>\n");

  std::ifstream infile;
  // Throw exception on errors.
  infile.exceptions (std::ios::failbit | std::ios::badbit | std::ios::eofbit);
  // Open input file.
  infile.open (argv[1]);
  // Read the graph from the input file.
  read_graph (infile);
  
  // Print the graph.
  if (rank == 0)
    {
      std::cout << std::endl;
      for (unsigned i = 1; i <= N; ++i)
        {
          for (unsigned j = 1; j <= N; ++j)
            std::cout << (*graph) (j, i) << " ";
          std::cout << std::endl;
        }
    }

  initialize ();
  // Divide the work and send it to the other nodes.
  if (rank == 0)
    {
      unsigned const amount = edge_count / worldsize;
      for (int i = 1; i < worldsize; ++i)
        {
          std::list<StackElem> tmplist;
          StackElem & el = stack.back ();
          unsigned cnt = 0;
          while (!(cnt > amount || el.next >= max_edges))
            {
              // Save the number of the soon to be used edge.
              unsigned const inserted = el.next;
              // Adjust the pointer to next edge.
              el.next = edgeinfo[inserted].next;

              /* el:    [1 0 0 ... 0]
                 |
                 v
                 newel: [1 1 0 ... 0]  */
              StackElem newel (el);
              // Increment used edges counter for the new element.
              newel.used_edges = el.used_edges + 1;
      
              /*
                Check whether the new state is a valid one, i.e. that adding
                the edge with either colour won't result in a monochromatic
                triangle.  */
              std::pair<bool, bool> const res = check_triangle (newel, 
                                                                inserted);
              assert (!res.first && !res.second);

              /*
                Push newel with one colour onto DFS stack if it doesn't create
                a monochromatic triangle.  */
              if (res.first == false)
                {
                  newel.colour[inserted] = false;
                  newel.valid[inserted] = true;
                  tmplist.push_back (newel);
                  states += 1;
                  cnt += 1;
                }

              /*
                Push newel with the other colour onto DFS stack if it doesn't
                create monochromatic triangle.  */
              if (res.second == false)
                {
                  newel.colour[inserted] = true;
                  newel.valid[inserted] = true;
                  tmplist.push_back (newel);
                  states += 1;
                  cnt += 1;
                }
            }
          log << "generated " << cnt << " elements" << std::endl;
          // Deny the request if we were not able to generate valid states.
          assert (cnt);
  
          /* Send the amount to requester. */
          mpi_obuffer obuf;
          // Message type and number of elements.
          obuf << TYPE_STKELEM << cnt;
          for (std::list<StackElem>::reverse_iterator it = tmplist.rbegin ();
               it != tmplist.rend (); ++it)
            {
              serialize_stackelem (obuf, *it);
            }
          // Finally the send.
          CHECKMPI (MPI_Send (obuf.data (), obuf.get_pos (), MPI_PACKED, i, 
                              TAG_NEEDS_ATTENTION, MPI_COMM_WORLD), 
                    "MPI_Send");
        }
    }
  /* Synchronize before start of the computation. */
  MPI_Barrier (MPI_COMM_WORLD);
  //unsigned old100k = 0;
  t1 = MPI_Wtime ();
  while (1)
    {
      char dnr;
      int flag = 0;

      /* Probe for incoming messages and process them. */
      while (1)
        {
          /* Probe for and process messages only after we have done some work
             or the stack is empty. */
          if (batch_size >= BATCHSIZE || stack.empty ())
            batch_size = 0;
          else
            break;

          flag = 0;
          CHECKMPI (MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
                                &flag, &status), "MPI_Iprobe");
          if (flag)
            {
              CHECKMPI (MPI_Recv (recv_buf, recv_buf_len, MPI_PACKED,
                                  MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
                                  &status), "MPI_Recv");
              mpi_ibuffer ibuf (recv_buf, recv_buf_len);
              switch (status.MPI_TAG)
                {
                case TAG_NEEDS_ATTENTION:
                  process_priority_message (ibuf);
                  break;
                default:
                  error ("Unknown TAG_*!!!");
                }
              continue;
            }
          else
            break;
        }

      /* Are we out of work? */
      if (stack.empty ())
        {
          log << "out of work" << std::endl;
          /* Deny any requests for work. */
          wouldgive = false;

          /* First do the right thing with tokens. */
          do_tokens ();

          if (rank != 0)
            /* Request a donor from 0. */
            dnr = request_donor ();
          else
            {
              dnr = donor;
              donor = (donor + 1) % worldsize;
            }
          if (dnr != rank)
            /* Send request to obtained donor and read results. */
            request_work (dnr);

          /* At this point we should have some work to give. */
          wouldgive = true;

          continue;
        }

      // Select top element of DFS stack.
      StackElem & el = stack.back ();

      // Is this the solution?
      if (el.used_edges == edge_count)
        {
          result_element = new StackElem (el);
          end_computation (true);
        }
        
      // Is it possible to go deeper in DFS tree?
      if (el.next >= max_edges)
        {
          // No, return instead.
          stack.pop_back ();
          continue;
        }
      
      assert (edges[el.next]);
      
      // Save the number of the soon to be used edge.
      unsigned const inserted = el.next;
      // Adjust the pointer to next edge.
      el.next = edgeinfo[inserted].next;

      /* el:    [1 0 0 ... 0]
         |
         v
         newel: [1 1 0 ... 0]  */
      StackElem newel (el);
      // Increment used edges counter for the new element.
      newel.used_edges = el.used_edges + 1;
      
      /*
        Check whether the new state is a valid one, i.e. that adding the edge
        with either colour won't result in a monochromatic triangle.  */
      std::pair<bool, bool> const res = check_triangle (newel, inserted);

      /*
        Push newel with one colour onto DFS stack if it doesn't create a
        monochromatic triangle.  */
      if (res.first == false)
        {
          newel.colour[inserted] = false;
          newel.valid[inserted] = true;
          stack.push_back (newel);
          states += 1;
          batch_size += 1;
        }

      /*
        Push newel with the other colour onto DFS stack if it doesn't create
        monochromatic triangle.  */
      if (res.second == false)
        {
          newel.colour[inserted] = true;
          newel.valid[inserted] = true;
          stack.push_back (newel);
          states += 1;
          batch_size += 1;
        }
    }

  // No solution has been found after a search through the whole state space.
  std::cout << "No solution has been found." << std::endl;

  return 0;
}
catch (std::exception & e)
{
  std::cerr << e.what () << std::endl;
}
catch (...)
{
  error ("Unknown exception has been caught.");
}
