/*******************************************************************************
  author       : Wasinee R.
  created      : DATE: 07.04.2003
  $Id$

__copyright__

*******************************************************************************/

typedef enum eDFSCOLORS {GRAY=0, BLACK=1, WHITE=2, NONE=-1} DFSCOLORS;

#include "mes.h"
#include "ghmm.h"
#include "model.h"

typedef struct local_store_topo {
  int     *topo_order;
  int     topo_order_length;
  int     *queue;
  int     head, tail;
} local_store_topo;
 
static local_store_topo *topo_alloc(model *mo, int len);
static int topo_free(local_store_topo **v, int n, int len);


/*----------------------------------------------------------------------------*/
static local_store_topo *topo_alloc(model *mo, int len) {
#define CUR_PROC "sdtopo_alloc"
  local_store_topo* v = NULL;
  int j;

  if (!m_calloc(v, 1)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->queue, mo->N)) {mes_proc(); goto STOP;}

  v->topo_order_length = 0;
  v->head =0; /* initialize static queue (array implementation) */
  v->tail =0;
  if (!m_calloc(v->topo_order, mo->N)) {mes_proc(); goto STOP;}

  return(v);
STOP:
  topo_free(&v, mo->N, len);
  return(NULL);
#undef CUR_PROC
} 


/*----------------------------------------------------------------------------*/
static int topo_free(local_store_topo **v, int n, int len) {
#define CUR_PROC "sdviterbi_free"
  int j;
  mes_check_ptr(v, return(-1));
  if( !*v ) return(0);
  (*v)->head=0; (*v)->tail=0;
  m_free((*v)->topo_order);
  m_free((*v)->queue);
  m_free(*v);
  return(0);
#undef CUR_PROC
} 

/*---------------------------------------------------------------------------*/
/*
 * Implementation of DFSVisit with recursion (WS)
 *
 */
static void model_DFSVisit(model *c_model, int nodev, int *timevisit, int *parents, int *colors, int **edge_classes)
{
#define CUR_PROC "DFSVisit"
  int i,w;
  colors[nodev]=GRAY;
  ++(*timevisit);
  for(i=0; i < c_model->s[nodev].out_states; i++) { 
    w=c_model->s[nodev].out_id[i]; // Explore edge (v,w)
    if ( edge_classes[nodev][w] == NONE ) { // First exploration
      edge_classes[nodev][w] = colors[w];
      /*fprintf(stderr, " %d edge (%s, %s)\n", (int) colors[w], c_model->s[nodev].label, c_model->s[w].label); */
    }
    if ( colors[w] == WHITE ) {
      parents[w]=nodev;
      model_DFSVisit(c_model, w, timevisit, parents, colors, edge_classes);
    } 
  }
  colors[nodev]=BLACK; // finished  
  ++(*timevisit);
#undef CUR_PROC
}


/** Return classification of edges in the model graph 
    WHITE EDGE => edges in the DFS search tree
    GRAY EDGE  => edges that form cycles which must be removed 
                  before running topological sort
*/
int** model_DFS(model *c_model)
{
#define CUR_PROC "model_DFS"
  int i,j,initials;
  int timevisit=0;
  int *colors;
  int *parents;
  int **edge_classes; 

  colors=(int*)calloc(c_model->N,sizeof(int));
  parents=(int*)calloc(c_model->N,sizeof(int));
  edge_classes=(int**)calloc(c_model->N,sizeof(int*));
  for(i=0; i < c_model->N; i++) {
    edge_classes[i]=(int*)calloc(c_model->N,sizeof(int));
  }
  for(i=0; i < c_model->N; i++) {
    if ( c_model->s[i].pi == 1.0) initials=i; /* assuming only one initial state */
    colors[i]=WHITE;
    parents[i]=-1;
    for(j=0; j < c_model->N; j++) {
      edge_classes[i][j] = NONE;
    }
  }
  
  model_DFSVisit(c_model, initials, &timevisit,parents, colors, edge_classes);
  for(i=0; i < c_model->N; i++) { /* handle subtrees */
    if ( colors[i] == WHITE ) {
      model_DFSVisit(c_model, i, &timevisit,parents, colors, edge_classes);
    }
  }

  m_free(colors);
  m_free(parents);
  return edge_classes;
STOP:
  return NULL;
#undef CUR_PROC 
}

static void __topological_sort( model *c_model, local_store_topo *v, int** edge_classes)
{
  int i,j;
  int nodeu, nodew, dels_cnt;
  int *indegrees = (int*) calloc(c_model->N,sizeof(int));
 
  for(i=0; i < c_model->N; i++) {
    indegrees[i]=c_model->s[i].in_states;
  }
  for(i=0; i < c_model->N; i++) { /* don't consider back edges in top sort */
    for(j=0; j < c_model->N; j++) {
      if (edge_classes[i][j] == GRAY) {
        indegrees[j]--;
        /* fprintf(stderr, "BACK edge (%s, %s)\n", c_model->s[i].label, c_model->s[j].label); */
      }
    }
  }
  /* Create a queue q and enqueue all nodes of indegree 0. */
  v->head=0;  v->tail=0;
  for(i=0; i < c_model->N; i++) {
    if ( indegrees[i] == 0 ) v->queue[v->tail++]=i;
  }
  dels_cnt=0;
  while( v->head != v->tail ) {
    nodeu=v->queue[v->head++]; /* dequeue */
    if (c_model->silent[nodeu]) {
      v->topo_order[dels_cnt++]=nodeu; /* append it to the list */
      /* printf("Silent state: %s\n", c_model->s[nodeu].label); */
    }
    for(i=0; i < c_model->s[nodeu].out_states; i++) {
      nodew=c_model->s[nodeu].out_id[i];
      if ( edge_classes[nodeu][nodew] != GRAY ) {
        indegrees[nodew]--;
        if (nodew != nodeu && indegrees[nodew]==0) {
          v->queue[v->tail++]=nodew; /* enqueue */
        }
      }
    }
  }
  v->topo_order_length=dels_cnt;
  free(indegrees); 
}

/*----------------------------------------------------------------------------*/
void model_topo_ordering(model *mo)
{
#define CUR_PROC "model_topo_ordering"
  int i;
  local_store_topo *v;
  int **edge_cls;

  v = topo_alloc(mo, 1);
  if (!v) { mes_proc(); goto STOP; }
 
  edge_cls=model_DFS( mo);
  __topological_sort( mo, v, edge_cls );
  mo->topo_order_length=v->topo_order_length;
  if (!m_calloc(mo->topo_order, mo->topo_order_length)) {mes_proc(); goto STOP;}
 
  for(i=0; i < v->topo_order_length; i++) {
    mo->topo_order[i] = v->topo_order[i];
  }
  
  /* fprintf(stderr,"Ordering silent states....\n\t");
  for(i=0; i < mo->topo_order_length; i++) {
    fprintf(stderr, "%d, ", mo->topo_order[i]);
  }
  fprintf(stderr,"\n\n");  */
  
  matrix_i_free(&edge_cls, mo->N);
  topo_free(&v, mo->N, 1);
STOP:
 i = 0;
 #undef CUR_PROC
}