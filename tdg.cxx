#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <utility>
#include <stack>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cerrno>
#include <cassert>
#include "matrix.hxx"


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

  StackElem (unsigned width, unsigned nxt, unsigned used)
    : colour (width), valid (width), next (nxt), used_edges (used)
  { }
};


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


//! Number of vertices in graph.
unsigned N = 0;
//! Maximal number of edges in graph.
unsigned max_edges = 0;
//! Number of edges in graph.
unsigned edge_count = 0;
//! First edge.
unsigned first_edge = 0;
//! Upper triangle matrix representing the graph.
Graph * graph = 0;

/*! 
  Vector of sets of neighbouring vertices for each vertex. For fast triangle
  testing.  */
Neighbours neighbours;

/*!
  \brief Bitmap mapping integers for all possible edges to booleans,
  edges[num] == true iff the edge is present in the graph.  */
ValidEdges edges;

/*!
  \brief Map from edge numbers to pairs of vertices. edgeinfo[num] is only
  valid iff edges[num] == true.  */
std::vector<EdgeInfo> edgeinfo;
//! DFS algoritm's stack.
std::stack<StackElem> stack;


/*!
  \brief Prints msg and possible message for errno to stderr and exits 
  with EXIT_FAILURE.
  \param msg User supplied message.  */
void 
error (const char * msg)
{
  if (! errno)
    std::cerr << msg << std::endl;
  else
    std::cerr << msg << ": " << strerror (errno) << std::endl;
  exit (EXIT_FAILURE);
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
  assert (el.valid[edge]);

  unsigned const v1 = edgeinfo[edge].v1;
  unsigned const v2 = edgeinfo[edge].v2;
  // Loops cannot create triangles with other edges.
  if (v1 == v2)
    return std::make_pair (false, false);
  
  /*
    Intersection of sets of neighbours of vertices that make edge gives us
    candidates vertices for triangle.  */
  //typedef std::list<unsigned> CommonVertices;
  typedef std::vector<unsigned> CommonVertices;
  CommonVertices common;
  common.reserve (N - 1);
  std::set_intersection (neighbours[v1-1].begin (), neighbours[v1-1].end (),
                         neighbours[v2-1].begin (), neighbours[v2-1].end (),
                         std::back_inserter (common));

  /*
    Test if candidate edges have the same colour. Loop until we either test
    all vertices or until we determine that adding vertex of either colour
    would create a triangle.  */
  std::pair<bool, bool> res (false, false);
  for (CommonVertices::const_iterator it = common.begin ();
       it != common.end () || (res.first && res.second);
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
   \brief Read graph from file and initializes all necessary variable.  
   
   \param infile Input stream.  */
void
read_graph (std::istream & infile)
{
  // Read graph's size from the input file.
  infile >> N;
  // Allocate and initialize structures.
  init_structures ();
  // Read the graph from input the file.
  unsigned last = 0;
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= N; ++j)
      {
        unsigned val;
        // Read edge from file.
        infile >> val;

        if (i > j)
          continue;

        // Update graph.
        (*graph)(j, i) = val;
        // Is there an edge here?
        if (val)
          {
            // Update number of edges.
            edge_count += 1;

            // Record the neighbourness of the two vertices.
            if (j != i)
              {
                neighbours[j-1].insert (i);
                neighbours[i-1].insert (j);
              }
            
            unsigned const idx = index_for_xy (j, i);

            // Is this the first edge that we see?
            if (!last)
              {
                // Yes. No successor edge yet, next loop will set it.
                edgeinfo[idx] = EdgeInfo (j, i);
                first_edge = idx;
              }
            else
              {
                // Set successor of the previous edge.
                edgeinfo[last].next = idx;
                // Successor will be set successor edge.
                edgeinfo[idx] = EdgeInfo (j, i);

                // Debug.
                std::cerr << " successor edge " << idx << std::endl;
              }
            // Update bitmap of valid edges.
            edges[idx] = true;
            // Record this edge for the next loop.
            last = idx;

            // Debug.
            std::cerr << "edge " << idx << " vertices (" << j << ", " << i << ")";
          }
      }
  // The last edge doesn't have any successor.
  edgeinfo[last].next = std::numeric_limits<unsigned>::max ();
  // We have counted each edge twice. Fixing it now.
  //edge_count /= 2;
}


/*!
   \brief Prints the solution represented by stack element.  

   \param solution Stack element representing the solution.  */
void
print_solution (StackElem const & solution)
{
  std::cout << "Solution has been found:" << std::endl;
  
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
  std::cout << "N: " << N << std::endl;
  std::cout << "edges: " << edge_count << std::endl;

  exit (EXIT_SUCCESS);
}


int
main (int argc, char * argv[])
try
{
  // Some basic checks and initialization.
  if (argc != 2)
    error ("Syntax: tdg <input graph>");

  std::ifstream infile;
  // Throw exception on errors.
  infile.exceptions (std::ios::failbit | std::ios::badbit | std::ios::eofbit);
  // Open input file.
  infile.open (argv[1]);
  // Read the graph from the input file.
  read_graph (infile);
  
  // Print the graph.
  std::cout << std::endl;
  for (unsigned i = 1; i <= N; ++i)
    {
      for (unsigned j = 1; j <= N; ++j)
        std::cout << (*graph) (j, i) << " ";
      std::cout << std::endl;
    }

  // Put initial state on the stack.
  stack.push (StackElem (max_edges, first_edge, 0));
  // Main loop.
  while (! stack.empty ())
    {
      // Select top element of DFS stack.
      StackElem & el = stack.top ();

      // Is this the solution?
      if (el.used_edges == edge_count)
        print_solution (el);
        
      // Is it possible to go deeper in DFS tree?
      if (el.next >= max_edges)
        {
          // No, return instead.
          stack.pop ();
          continue;
        }
      
      assert (edges[el.next]);
      
      // Save the number of the soon to be used edge.
      unsigned const inserted = el.next;
      // Adjust the pointer to next edge.
      el.next = edgeinfo[el.next].next;

      /* el:    [1 0 0 ... 0]
         |
         v
         newel: [1 1 0 ... 0]  */
      StackElem newel (el);
      // Add inserted edge to the set of used edges.
      newel.valid[inserted] = true;
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
          stack.push (newel);
        }

      /*
        Push newel with the other colour onto DFS stack if it doesn't create
        monochromatic triangle.  */
      if (res.second == false)
        {
          newel.colour[inserted] = true;
          stack.push (newel);
        }
    }

  // No solution has been found after a search through the whole state space.
  std::cout << "No solution has been found." << std::endl;

  return 0;
}
catch (std::exception & e)
{
  error (e.what ());
}
catch (...)
{
  error ("Unknown exception has been caught.");
}
