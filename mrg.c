#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "config.h"
#include "bitmap.h"
#include "matrix.h"
#include "list.h"


struct _stkelem_t 
{
  /* Weight of this cut. */
  int weight;
  /* Offset of the rightmost 1 in the future new element
     generated from this element.*/
  unsigned next;
  /* Representation of X and Y sets. */
  bitmap_t * set;
};
typedef struct _stkelem_t stkelem_t;


/* Number of nodes. */
unsigned N = 0;
/* Stack for DFS algorithm. */
list_t * stack;
/* */
trimatrix_t * graph;
/* Matrix of edges' weights. */
wtrimatrix_t * weights;
/* Best solution. */
stkelem_t * best;


/**
   Initializes new DFS stack element.
   @param se pointer to stack element
   @param width width of bitmap/set
   @param weight weight of cut in this step
   @param rightmost offset of the rightmost 1 in bitmap/set
   @param last offset of the last 1 in the last generated element
*/
inline
stkelem_t *
stkelem_init (stkelem_t * se, unsigned width, unsigned weight, 
              unsigned next)
{
  if (width == 0 || next > width - 1)
    abort ();
  
  se->set = bitmap_new (width);
  if (! se->set)
    {
      free (se);
      return NULL;
    }
  se->weight = weight;
  se->next = next;
  return se;
}


/**
   Allocates and initializes new DFS stack element.
   @param width width of bitmap/set
   @param weight weight of cut in this step
   @param rightmost offset of the rightmost 1 in bitmap/set
   @param next offset of the rightmost 1 of future generated element
*/
stkelem_t * 
stkelem_new (unsigned width, unsigned weight, unsigned next)
{
  stkelem_t * se;

  se = malloc (sizeof (stkelem_t));
  if (! se)
    return NULL;
  if (! stkelem_init (se, width, weight, next))
    {
      free (se);
      return NULL;
    }
  return se;
}


/**
   Clones DFS stack element.
   @param se stack element
   @return copy of the stack element
*/
stkelem_t * 
stkelem_clone (const stkelem_t * se)
{
  stkelem_t * newse;
  
  newse = malloc (sizeof (stkelem_t));
  if (! newse)
    return NULL;
  newse->set = bitmap_clone (se->set);
  if (! newse->set)
    {
      free (newse);
      return NULL;
    }
  newse->weight = se->weight;
  newse->next = se->next;
  return newse;
}


/**

*/
inline
void 
stkelem_destroy (const stkelem_t * se)
{
  bitmap_delete (se->set);
}


/**
 
*/
void
stkelem_delete (stkelem_t * se)
{
  stkelem_destroy (se);
  free (se);
}


/**
   Prints msg and possible message for errno to stderr and 
   exits with EXIT_FAILURE.
   @param msg user supplied message
*/
void 
error (const char * msg)
{
  if (! errno)
    fprintf (stderr, "%s\n", msg);
  else
    fprintf (stderr, "%s: %s\n", msg, strerror (errno));
  exit (EXIT_FAILURE);
}


/**
   Initializes stack for DFS algoritm.
*/
void 
initialize_stack (void)
{
  stkelem_t * el = stkelem_new (N, 0, 0);
  
  if (! el)
    error ("Memory allocation failure");
  if (! list_pushback (stack, el))
    error ("list_pushback()");
}


/**
   Global initialization of computation.
*/
void
initialize (void)
{
  initialize_stack ();
  best = stkelem_new (N, INT_MAX, 0);
  if (! best)
    error ("Memory allocation failure");
}


/**
   Updates weight of cut when we move one node from set X to Y.
   @param el element of DFS tree to update
   @param node node that has been moved from X to Y
   @return true if the cut has weight 1, false otherwise
*/
int 
update_weight (stkelem_t * el, unsigned node)
{
  unsigned i;

  if (node == 0)
    abort ();

  for (i = 1; i <= N; ++i)
    {
      if (i == node)
        continue;
      if (trimatrix_get (graph, node, i))
        {
          /* Is node i in set Y? */
          if (bitmap_getbit (el->set, i-1))
            /* Substract weight of edges whose end nodes are now
               both in Y from the weight of the cut. */
            el->weight -= wtrimatrix_get (weights, node, i);
          else
            /* Add weight of edges whose end nodes are now one in
               the set X and the other in the set Y. */
            el->weight += wtrimatrix_get (weights, node, i);
        }
    }
  
  if (el->weight < best->weight && el->weight > 0)
    best = el;

  if (el->weight == 1)
    return 1;
  else
    return 0;
}


/**
   Generates next level of DFS tree from element el and pushes it 
   onto DFS stack.
   @param el element
   @return true if the next element was successfully generated,
   false otherwise.
*/
int
generate_depth (stkelem_t * el)
{
  stkelem_t * newel;

  /* Is it possible to go deeper in DFS tree? */
  if (el->next < N)
    {
      /* el:    [1 0 0 ... 0]
         |        
         v        
         newel: [1 1 0 ... 0] */
      newel = stkelem_clone (el);
      if (! newel)
        error ("Memory allocation failure");
      bitmap_setbit (newel->set, el->next);
      newel->next = el->next + 1;
      el->next += 1;
      /* Push newel onto DFS stack. */
      if (! list_push (stack, newel))
        error ("list_push()");
      return 1;
    }
  else
    return 0;
}
  

int 
main (int argc, char * argv[])
{
  int ret;
  unsigned i, j;
  FILE * infile;

  /* Some basic checks and initialization. */
  if (argc != 2)
    error ("Syntax: mrg <input graph>");
  srandom (time (NULL));
  /* Open input file and read graph's dimension. */
  infile = fopen (argv[1], "r");
  if (! infile)
    error ("fopen()");
  ret = fscanf (infile, "%u", &N);
  if (ret < 1)
    error ("fscanf()");
  /* Allocate structures. */
  stack = list_new ();
  graph = trimatrix_new (N);
  weights = wtrimatrix_new (N);
  if (! stack || ! graph || ! weights)
    error ("Memory allocation failure");
  /* Read graph from file. */
  for (i = 1; i <= N; ++i)
    for (j = 1; j <= N; ++j)
      {
        unsigned val;
        ret = fscanf (infile, "%u", &val);
        if (ret < 1)
          error ("fscanf()");
        trimatrix_set (graph, i, j, val);
        if (val)
          wtrimatrix_set (weights, i, j, random () % 255 + 1);
      }

  /* Do the actual work here.  */
  initialize ();
  while (1)
    {
      listelem_t * it;
      stkelem_t * el;
      
      el = list_first (stack, &it);
      if (! el)
        break;
      /* Move deeper in DFS tree if possible. */
      if (generate_depth (el))
        {
          /* Get the newly generated element. */
          el = list_first (stack, &it);
          /* Update weight of a new cut.
             Is this a cut of weight 1? 
             Note: el->next because nodes are numbered from 1. */
          if (update_weight (el, el->next)) 
            /* It is, we are done. */
            break;
          else
            continue;
        }
      else
        {
          stkelem_t * se = list_pop (stack);
          if (se != best)
            stkelem_delete (se);
        }
    }

  /* Print out the solution. */
  printf ("Weight of the best solution: %d\n", best->weight);
  printf ("Set X:");
  for (i = 0; i < bitmap_size (best->set); ++i)
    {
      int b = bitmap_getbit (best->set, i);
      if (! b)
	printf (" %d", i+1);
    }
  printf ("\n");
  printf ("Set Y:");
  for (i = 0; i < bitmap_size (best->set); ++i)
    {
      int b = bitmap_getbit (best->set, i);
      if (b)
	printf (" %d", i+1);
    }
  printf ("\n");
  
  exit (EXIT_SUCCESS);
}
