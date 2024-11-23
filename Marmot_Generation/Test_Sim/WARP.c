/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: WARP.instr (template_simple)
 * Date:       Mon Nov 04 14:54:42 2024
 * File:       WARP.c
 * CFLAGS=
 */

#define MCCODE_STRING "McStas 3.4 - Sep. 19, 2023"
#define FLAVOR        "mcstas"
#define FLAVOR_UPPER  "MCSTAS"

#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED

#include <string.h>

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];
#define MCCODE_BASE_TYPES

#ifndef MC_NUSERVAR
#define MC_NUSERVAR 10
#endif

/* Particle JUMP control logic */
struct particle_logic_struct {
int dummy;
};

struct _struct_particle {
  double x,y,z; /* position [m] */
  double vx,vy,vz; /* velocity [m/s] */
  double sx,sy,sz; /* spin [0-1] */
  int mcgravitation; /* gravity-state */
  void *mcMagnet;    /* precession-state */
  int allow_backprop; /* allow backprop */
  /* Generic Temporaries: */
  /* May be used internally by components e.g. for special */
  /* return-values from functions used in trace, thusreturned via */
  /* particle struct. (Example: Wolter Conics from McStas, silicon slabs.) */
  double _mctmp_a; /* temp a */
  double _mctmp_b; /* temp b */
  double _mctmp_c; /* temp c */
  unsigned long randstate[7];
  double t, p;    /* time, event weight */
  long long _uid;  /* event ID */
  long _index;     /* component index where to send this event */
  long _absorbed;  /* flag set to TRUE when this event is to be removed/ignored */
  long _scattered; /* flag set to TRUE when this event has interacted with the last component instance */
  long _restore;   /* set to true if neutron event must be restored */
  long flag_nocoordschange;   /* set to true if particle is jumping */
  struct particle_logic_struct _logic;
};
typedef struct _struct_particle _class_particle;

_class_particle _particle_global_randnbuse_var;
_class_particle* _particle = &_particle_global_randnbuse_var;

#pragma acc routine
_class_particle mcgenstate(void);
#pragma acc routine
_class_particle mcsetstate(double x, double y, double z, double vx, double vy, double vz,
			   double t, double sx, double sy, double sz, double p, int mcgravitation, void *mcMagnet, int mcallowbackprop);

extern int mcgravitation;      /* flag to enable gravitation */
#pragma acc declare create ( mcgravitation )
int mcallowbackprop;        
#pragma acc declare create ( mcallowbackprop )

_class_particle mcgenstate(void) {
  _class_particle particle = mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, mcgravitation, NULL, mcallowbackprop);
  return(particle);
}
/*Generated user variable handlers:*/

#pragma acc routine
double particle_getvar(_class_particle *p, char *name, int *suc);

#ifdef OPENACC
#pragma acc routine
int str_comp(char *str1, char *str2);
#endif

double particle_getvar(_class_particle *p, char *name, int *suc){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int s=1;
  double rval=0;
  if(!str_comp("x",name)){rval=p->x;s=0;}
  if(!str_comp("y",name)){rval=p->y;s=0;}
  if(!str_comp("z",name)){rval=p->z;s=0;}
  if(!str_comp("vx",name)){rval=p->vx;s=0;}
  if(!str_comp("vy",name)){rval=p->vy;s=0;}
  if(!str_comp("vz",name)){rval=p->vz;s=0;}
  if(!str_comp("sx",name)){rval=p->sx;s=0;}
  if(!str_comp("sy",name)){rval=p->sy;s=0;}
  if(!str_comp("sz",name)){rval=p->sz;s=0;}
  if(!str_comp("t",name)){rval=p->t;s=0;}
  if(!str_comp("p",name)){rval=p->p;s=0;}
  if(!str_comp("_mctmp_a",name)){rval=p->_mctmp_a;s=0;}
  if(!str_comp("_mctmp_b",name)){rval=p->_mctmp_b;s=0;}
  if(!str_comp("_mctmp_c",name)){rval=p->_mctmp_c;s=0;}
  if (suc!=0x0) {*suc=s;}
  return rval;
}

#pragma acc routine
void* particle_getvar_void(_class_particle *p, char *name, int *suc);

#ifdef OPENACC
#pragma acc routine
int str_comp(char *str1, char *str2);
#endif

void* particle_getvar_void(_class_particle *p, char *name, int *suc){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int s=1;
  void* rval=0;
  if(!str_comp("x",name)) {rval=(void*)&(p->x); s=0;}
  if(!str_comp("y",name)) {rval=(void*)&(p->y); s=0;}
  if(!str_comp("z",name)) {rval=(void*)&(p->z); s=0;}
  if(!str_comp("vx",name)){rval=(void*)&(p->vx);s=0;}
  if(!str_comp("vy",name)){rval=(void*)&(p->vy);s=0;}
  if(!str_comp("vz",name)){rval=(void*)&(p->vz);s=0;}
  if(!str_comp("sx",name)){rval=(void*)&(p->sx);s=0;}
  if(!str_comp("sy",name)){rval=(void*)&(p->sy);s=0;}
  if(!str_comp("sz",name)){rval=(void*)&(p->sz);s=0;}
  if(!str_comp("t",name)) {rval=(void*)&(p->t); s=0;}
  if(!str_comp("p",name)) {rval=(void*)&(p->p); s=0;}
  if (suc!=0x0) {*suc=s;}
  return rval;
}

#pragma acc routine
int particle_setvar_void(_class_particle *, char *, void*);

int particle_setvar_void(_class_particle *p, char *name, void* value){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int rval=1;
  if(!str_comp("x",name)) {memcpy(&(p->x),  value, sizeof(double)); rval=0;}
  if(!str_comp("y",name)) {memcpy(&(p->y),  value, sizeof(double)); rval=0;}
  if(!str_comp("z",name)) {memcpy(&(p->z),  value, sizeof(double)); rval=0;}
  if(!str_comp("vx",name)){memcpy(&(p->vx), value, sizeof(double)); rval=0;}
  if(!str_comp("vy",name)){memcpy(&(p->vy), value, sizeof(double)); rval=0;}
  if(!str_comp("vz",name)){memcpy(&(p->vz), value, sizeof(double)); rval=0;}
  if(!str_comp("sx",name)){memcpy(&(p->sx), value, sizeof(double)); rval=0;}
  if(!str_comp("sy",name)){memcpy(&(p->sy), value, sizeof(double)); rval=0;}
  if(!str_comp("sz",name)){memcpy(&(p->sz), value, sizeof(double)); rval=0;}
  if(!str_comp("p",name)) {memcpy(&(p->p),  value, sizeof(double)); rval=0;}
  if(!str_comp("t",name)) {memcpy(&(p->t),  value, sizeof(double)); rval=0;}
  return rval;
}

#pragma acc routine
int particle_setvar_void_array(_class_particle *, char *, void*, int);

int particle_setvar_void_array(_class_particle *p, char *name, void* value, int elements){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int rval=1;
  return rval;
}

#pragma acc routine
void particle_restore(_class_particle *p, _class_particle *p0);

void particle_restore(_class_particle *p, _class_particle *p0) {
  p->x  = p0->x;  p->y  = p0->y;  p->z  = p0->z;
  p->vx = p0->vx; p->vy = p0->vy; p->vz = p0->vz;
  p->sx = p0->sx; p->sy = p0->sy; p->sz = p0->sz;
  p->t = p0->t;  p->p  = p0->p;
  p->_absorbed=0; p->_restore=0;
}

#pragma acc routine
double particle_getuservar_byid(_class_particle *p, int id, int *suc){
  int s=1;
  double rval=0;
  switch(id){
  }
  if (suc!=0x0) {*suc=s;}
  return rval;
}

#pragma acc routine
void particle_uservar_init(_class_particle *p){
}

#define MC_EMBEDDED_RUNTIME
/* embedding file "mccode-r.h" */

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 3.4
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int numipar;
*   metadata_table_t metadata_table[];
*   int num_metadata;
*   char instrument_name[], instrument_source[];
*   int traceenabled, defaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <float.h>
#include <inttypes.h>
#include <stdint.h>
#ifdef OPENACC
#include <openacc.h>
#ifndef GCCOFFLOAD
#include <accelmath.h>
#else
#include <math.h>
#endif
#pragma acc routine
int noprintf();
#pragma acc routine
size_t str_len(const char *s);
#else
#include <math.h>
#endif

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#  define mcstatic
#else
#  define mcstatic
#endif

#ifdef __dest_os
#  if (__dest_os == __mac_os)
#    define MAC
#  endif
#endif

#ifdef __FreeBSD__
#  define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#  define NEED_STAT_H
#endif

#ifdef WIN32
#  define NEED_STAT_H
#  define NEED_TYPES_H
#endif

#ifdef NEED_STAT_H
#  include <sys/stat.h>
#endif

#ifdef NEED_TYPES_H
#  include <sys/types.h>
#endif

#ifndef MC_PATHSEP_C
#  ifdef WIN32
#    define MC_PATHSEP_C '\\'
#    define MC_PATHSEP_S "\\"
#  else  /* !WIN32 */
#    define MC_PATHSEP_C '/'
#    define MC_PATHSEP_S "/"
#  endif /* !WIN32 */
#endif /* MC_PATHSEP_C */

#ifdef WIN32
#  define mkdir(a,b) mkdir(a)
#  define getpid() NULL
#endif

/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#  define MCCODE_STRING "McStas 3.4 - Sep. 19, 2023"
#endif

#ifndef MCCODE_DATE
#  define MCCODE_DATE "Sep. 19, 2023"
#endif

#ifndef MCCODE_VERSION
#  define MCCODE_VERSION "3.4"
#endif

#ifndef MCCODE_NAME
#  define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#  define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_PARTICLE_CODE
#  define MCCODE_PARTICLE_CODE 2112
#endif

#ifndef MCCODE_LIBENV
#  define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#  define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#ifdef MAC
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#if (USE_MPI == 0)
#  undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#ifdef OPENACC  /* default is to disable signals with PGI/OpenACC */
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#ifndef OPENACC
#  ifndef USE_OFF  /* default is to enable OFF when not using PGI/OpenACC */
#    define USE_OFF
#  endif
#  ifndef CPUFUNNEL  /* allow to enable FUNNEL-mode on CPU */
#  ifdef FUNNEL      /* by default disable FUNNEL-mode when not using PGI/OpenACC */
#    undef FUNNEL
#  endif
#  endif
#endif

#if (NOSIGNALS == 0)
#  undef NOSIGNALS
#endif

/** Header information for metadata-r.c ----------------------------------------------------------------------------- */
struct metadata_table_struct { /* stores metadata strings from components */
  char * source;  // component name which provided the metadata
  char * name;    // the name of the metadata
  char * type;    // the MIME type of the metadata (free form, valid identifier)
  char * value;   // the metadata string contents
};
typedef struct metadata_table_struct metadata_table_t;
char * metadata_table_key_component(char* key);
char * metadata_table_key_literal(char * key);
int metadata_table_defined(int, metadata_table_t *, char *);
char * metadata_table_type(int, metadata_table_t *, char *);
char * metadata_table_literal(int, metadata_table_t *, char *);
void metadata_table_print_all_keys(int no, metadata_table_t * tab);
int metadata_table_print_all_components(int no, metadata_table_t * tab);
int metadata_table_print_component_keys(int no, metadata_table_t * tab, char * key);
/* -------------------------------------------------------------------------- Header information for metadata-r.c --- */

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_int,
    instr_type_string, instr_type_char,
    instr_type_vector, instr_type_double
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
  char *unit; /* expected unit for parameter; informational only */
};


#ifndef MCCODE_BASE_TYPES
typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];
#endif

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[];         /* list of instrument parameters */
extern int    numipar;                                    /* number of instrument parameters */
extern metadata_table_t metadata_table[];                 /* list of component-defined string metadata */
extern int    num_metadata;                               /* number of component-defined string metadata */
extern char   instrument_name[], instrument_source[]; /* instrument name and filename */
extern char  *instrument_exe;                           /* executable path = argv[0] or NULL */
extern char   instrument_code[];                        /* contains the initial 'instr' file */

#ifndef MC_ANCIENT_COMPATIBILITY
extern int traceenabled, defaultmain;
#endif
#endif


/* Useful macros ============================================================ */


/* SECTION: Dynamic Arrays */
typedef int* IArray1d;
IArray1d create_iarr1d(int n);
void destroy_iarr1d(IArray1d a);

typedef int** IArray2d;
IArray2d create_iarr2d(int nx, int ny);
void destroy_iarr2d(IArray2d a);

typedef int*** IArray3d;
IArray3d create_iarr3d(int nx, int ny, int nz);
void destroy_iarr3d(IArray3d a);

typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);


/* MPI stuff */
#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#ifndef OPENACC    /* ... but only if we are not also running on GPU */
#undef NOSIGNALS
#endif
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */


#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
#pragma acc routine
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount-1 */

/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
char  *mcsig_message;
#define SIG_MESSAGE(msg) mcsig_message=(char *)(msg);
#else
#define SIG_MESSAGE(...)
#endif /* !NOSIGNALS */


/* Useful macros and constants ============================================== */


#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif


#  ifndef M_E
#    define M_E        2.71828182845904523536  // e
#  endif
#  ifndef M_LOG2E
#    define M_LOG2E    1.44269504088896340736  //  log2(e)
#  endif
#  ifndef M_LOG10E
#    define M_LOG10E   0.434294481903251827651 //  log10(e)
#  endif
#  ifndef M_LN2
#    define M_LN2      0.693147180559945309417 //  ln(2)
#  endif
#  ifndef M_LN10
#    define M_LN10     2.30258509299404568402  //  ln(10)
#  endif
#  ifndef M_PI
#    define M_PI       3.14159265358979323846  //  pi
#  endif
#  ifndef PI
#    define PI       M_PI                      //  pi - also used in some places
#  endif
#  ifndef M_PI_2
#    define M_PI_2     1.57079632679489661923  //  pi/2
#  endif
#  ifndef M_PI_4
#    define M_PI_4     0.785398163397448309616 //  pi/4
#  endif
#  ifndef M_1_PI
#    define M_1_PI     0.318309886183790671538 //  1/pi
#  endif
#  ifndef M_2_PI
#    define M_2_PI     0.636619772367581343076 //  2/pi
#  endif
#  ifndef M_2_SQRTPI
#    define M_2_SQRTPI 1.12837916709551257390  //  2/sqrt(pi)
#  endif
#  ifndef M_SQRT2
#    define M_SQRT2    1.41421356237309504880  //  sqrt(2)
#  endif
#  ifndef M_SQRT1_2
#    define M_SQRT1_2  0.707106781186547524401 //  1/sqrt(2)
#  endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


#define UNSET nan("0x6E6F74736574")
int nans_match(double, double);
int is_unset(double);
int is_valid(double);
int is_set(double);
int all_unset(int n, ...);
int all_set(int n, ...);
int any_unset(int n, ...);
int any_set(int n, ...);


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) (instrument->_position_absolute[index])
#define POS_R_COMP_INDEX(index) (instrument->_position_relative[index])

/* setting parameters based COMP_GETPAR (returned as pointer)         */
/* compname must be given as a string, type and par are symbols.      */
#define COMP_GETPAR3(type, compname, par) \
    &( ((_class_ ## type ##_parameters *) _getvar_parameters(compname))->par )
/* the body of this function depends on component instances, and is cogen'd */
void* _getvar_parameters(char* compname);

int _getcomp_index(char* compname);

/* Note: The two-stage approach to COMP_GETPAR is NOT redundant; without it,
* after #define C sample, COMP_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of COMP_GETPAR requires that we use sometimes bare names...
* NOTE: This can ONLY be used in instrument descriptions, not components.
*/
#define COMP_GETPAR2(comp, par) (_ ## comp ## _var._parameters.par)
#define COMP_GETPAR(comp, par) COMP_GETPAR2(comp,par)

#define INSTRUMENT_GETPAR(par) (_instrument_var._parameters.par)

/* Current component name, index, position and orientation */
/* These macros work because, using class-based functions, "comp" is usually
*  the local variable of the active/current component. */
#define INDEX_CURRENT_COMP (_comp->_index)
#define NAME_CURRENT_COMP (_comp->_name)
#define TYPE_CURRENT_COMP (_comp->_type)
#define POS_A_CURRENT_COMP (_comp->_position_absolute)
#define POS_R_CURRENT_COMP (_comp->_position_relative)
#define ROT_A_CURRENT_COMP (_comp->_rotation_absolute)
#define ROT_R_CURRENT_COMP (_comp->_rotation_relative)

#define NAME_INSTRUMENT (instrument->_name)


/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define DEBUG_INSTR() if(!mcdotrace); else { printf("INSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", instrument_name, instrument_source); }
#define DEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  printf("Component %30s AT (%g,%g,%g)\n", name, c.x, c.y, c.z); \
  }
#define DEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define DEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define DEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define DEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define DEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define DEBUG_INSTR()
#define DEBUG_COMPONENT(name,c,t)
#define DEBUG_INSTR_END()
#define DEBUG_ENTER()
#define DEBUG_COMP(c)
#define DEBUG_LEAVE()
#define DEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz);
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz);
void mcdis_sphere(double x, double y, double z, double r, int N);


/* random number generation. ================================================ */

/* available random number generators */
#define _RNG_ALG_MT         1
#define _RNG_ALG_KISS       2

/* selection of random number generator */
#ifndef RNG_ALG
#  define RNG_ALG  _RNG_ALG_KISS
#endif


#if RNG_ALG == _RNG_ALG_MT  // MT (currently not functional for gpu)
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define randstate_t unsigned long // this could be anything
#  define RANDSTATE_LEN 1
#  define srandom(seed) mt_srandom_empty()
#  define random() mt_random()
#  define _random() mt_random()
#elif RNG_ALG == _RNG_ALG_KISS  // KISS
#  ifndef ULONG_MAX
#    define ULONG_MAX ((unsigned long)0xffffffffffffffffUL)
#  endif
#  define MC_RAND_MAX ULONG_MAX
#  define randstate_t unsigned long
#  define RANDSTATE_LEN 7
#  define srandom(seed) kiss_srandom(_particle->randstate, seed)
#  define random() kiss_random(_particle->randstate)
#  define _random() kiss_random(state)
#endif

#pragma acc routine
double _randnorm2(randstate_t* state);


// component writers interface
#define randnorm() _randnorm2(_particle->randstate) // NOTE: can not use _randnorm on gpu
#define rand01() _rand01(_particle->randstate)
#define randpm1() _randpm1(_particle->randstate)
#define rand0max(p1) _rand0max(p1, _particle->randstate)
#define randminmax(p1, p2) _randminmax(p1, p2, _particle->randstate)
#define randtriangle() _randtriangle(_particle->randstate)

// Mersenne Twister rng
unsigned long mt_random(void);
void mt_srandom (unsigned long x);
void mt_srandom_empty();

// KISS rng
#pragma acc routine
unsigned long *kiss_srandom(unsigned long state[7], unsigned long seed);
#pragma acc routine
unsigned long kiss_random(unsigned long state[7]);

// Scrambler / hash function
#pragma acc routine seq
randstate_t _hash(randstate_t x);

// internal RNG (transforms) interface
#pragma acc routine
double _rand01(randstate_t* state);
#pragma acc routine
double _randpm1(randstate_t* state);
#pragma acc routine
double _rand0max(double max, randstate_t* state);
#pragma acc routine
double _randminmax(double min, double max, randstate_t* state);
#pragma acc routine
double _randtriangle(randstate_t* state);


#ifdef USE_OPENCL
#include "opencl-lib.h"
#include "opencl-lib.c"
#endif

#ifndef DANSE
int init(void);
int raytrace(_class_particle*);
int save(FILE *);
int finally(void);
int display(void);
#endif


/* GPU related algorithms =================================================== */

/*
*  Divide-and-conquer strategy for parallel sort absorbed last.
*/
#ifdef FUNNEL
long sort_absorb_last(_class_particle* particles, _class_particle* pbuffer, long len, long buffer_len, long flag_split, long* multiplier);
#endif
long sort_absorb_last_serial(_class_particle* particles, long len);


/* simple vector algebra ==================================================== */


#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
#pragma acc routine seq
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

#pragma acc routine seq
mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#pragma acc routine seq
mcstatic void norm_func(double *x, double *y, double *z);
#define NORM(x,y,z)	norm_func(&x, &y, &z)

#pragma acc routine seq
void normal_vec(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

#pragma acc routine
Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
#pragma acc routine
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
#pragma acc routine
Coords coords_add(Coords a, Coords b);
#pragma acc routine
Coords coords_sub(Coords a, Coords b);
#pragma acc routine
Coords coords_neg(Coords a);
#pragma acc routine
Coords coords_scale(Coords b, double scale);
#pragma acc routine
double coords_sp(Coords a, Coords b);
#pragma acc routine
Coords coords_xp(Coords b, Coords c);
#pragma acc routine
double coords_len(Coords a);
#pragma acc routine seq
void   coords_print(Coords a);
#pragma acc routine seq
mcstatic void coords_norm(Coords* c);

#pragma acc routine seq
void rot_set_rotation(Rotation t, double phx, double phy, double phz);
#pragma acc routine seq
int  rot_test_identity(Rotation t);
#pragma acc routine seq
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
#pragma acc routine seq
void rot_copy(Rotation dest, Rotation src);
#pragma acc routine seq
void rot_transpose(Rotation src, Rotation dst);
#pragma acc routine seq
Coords rot_apply(Rotation t, Coords a);

#pragma acc routine seq
void mccoordschange(Coords a, Rotation t, _class_particle *particle);
#pragma acc routine seq
void mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters
is no longer equal */

_class_particle mcgenstate(void);

// trajectory/shape intersection routines
#pragma acc routine seq
int inside_rectangle(double, double, double, double);
#pragma acc routine seq
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
      double vx, double vy, double vz, double dx, double dy, double dz);
#pragma acc routine seq
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
      double vx, double vy, double vz, double r, double h);
#pragma acc routine seq
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
      double vx, double vy, double vz, double r);
// second order equation roots
#pragma acc routine seq
int solve_2nd_order(double *t1, double *t2,
      double A,  double B,  double C);

// random vector generation to shape
// defines silently introducing _particle as the last argument
#define randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, radius) \
  _randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, radius, _particle)
#define randvec_target_rect_angular(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A) \
  _randvec_target_rect_angular(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, _particle)
#define randvec_target_rect_real(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, lx, ly, lz, order) \
  _randvec_target_rect_real(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, lx, ly, lz, order, _particle)
// defines forwarding to "inner" functions
#define randvec_target_sphere randvec_target_circle
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9) \
  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
// headers for randvec
#pragma acc routine seq
void _randvec_target_circle(double *xo, double *yo, double *zo,
  double *solid_angle, double xi, double yi, double zi, double radius,
  _class_particle* _particle);
#pragma acc routine seq
void _randvec_target_rect_angular(double *xo, double *yo, double *zo,
  double *solid_angle, double xi, double yi, double zi, double height,
  double width, Rotation A,
  _class_particle* _particle);
#pragma acc routine seq
void _randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
  double xi, double yi, double zi, double height, double width, Rotation A,
  double lx, double ly, double lz, int order,
  _class_particle* _particle);


// this is the main()
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif


/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *dirname             = NULL;      /* name of output directory */
static   char *siminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * siminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *siminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle);

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

/* embedding file "mcstas-r.h" */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char instrument_name[], instrument_source[];
*   int traceenabled, defaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  instrument.counter_AbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundent with mcstas.h */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER0 do {DEBUG_SCATTER(); SCATTERED++;} while(0)
#define SCATTER SCATTER0

#define JUMPTOCOMP(comp) mcneutron->_index = INDEX_COMP(comp);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  } while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(x, y, z, t, vx, vy, vz, \
               &sx, &sy, &sz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    x += vx*(dt); \
    y += vy*(dt); \
    z += vz*(dt); \
    t += (dt); \
    if (isnan(p) || isinf(p)) { ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { ABSORB; }\
    if (mcMagnet) /*printf("Spin precession gravity\n")*/; \
    x  += vx*(dt) + (Ax)*(dt)*(dt)/2; \
    y  += vy*(dt) + (Ay)*(dt)*(dt)/2; \
    z  += vz*(dt) + (Az)*(dt)*(dt)/2; \
    vx += (Ax)*(dt); \
    vy += (Ay)*(dt); \
    vz += (Az)*(dt); \
    t  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; ABSORB; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -vz, -z); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); z=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) { ABSORB; }; } \
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(vz == 0) { ABSORB; }; \
    mc_dt = -z/vz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { ABSORB; }; \
    mcPROP_DT(mc_dt); \
    z = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -vx, -x); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); x=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) { ABSORB; }; } \
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(vx == 0) { ABSORB; }; \
    mc_dt = -x/vx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { ABSORB; }; \
    mcPROP_DT(mc_dt); \
    x = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -vy, -y); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); y=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) { ABSORB; }; } \
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(vy == 0) { ABSORB; }; \
    mc_dt = -y/vy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { ABSORB; }; \
    mcPROP_DT(mc_dt); \
    y = 0; \
    DISALLOW_BACKPROP; \
  } while(0)


#ifdef DEBUG

#define DEBUG_STATE() if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define DEBUG_SCATTER() if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define DEBUG_STATE()
#define DEBUG_SCATTER()

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

/* embedding file "mccode-r.c" */

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#include <sys/stat.h>
#endif


#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int traceenabled = 0;
int defaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
#pragma acc declare create ( mcseed )
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
#pragma acc declare create ( mcdotrace )
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* OpenACC-related segmentation parameters: */
int vecsize = 128;
int numgangs = 7813;
long gpu_innerloop = 2147483647;

/* Monitor_nD list/buffer-size default */
long MONND_BUFSIZ = 1000000;

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
#ifdef MCDEFAULT_NCOUNT
mcstatic unsigned long long int mcncount             = MCDEFAULT_NCOUNT;
#else
mcstatic unsigned long long int mcncount             = 1000000;
#endif
#pragma acc declare create ( mcncount )
mcstatic unsigned long long int mcrun_num            = 0;
#pragma acc declare create ( mcrun_num )
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

#ifndef NX_COMPRESSION
#define NX_COMPRESSION NX_COMP_NONE
#endif

/* String nullification on GPU and other replacements */
#ifdef OPENACC
int noprintf() {
  return 0;
}

int str_comp(char *str1, char *str2) {
  while (*str1 && *str1 == *str2) {
    str1++;
    str2++;
  }
  return (*str1 - *str2);
}

size_t str_len(const char *s)
{
  size_t len = 0;
  if(s != NULL)
  {
    while(*s != '\0')
    {
      ++len;
      ++s;
    }
  }
  return len;
}

#endif

/* SECTION: Predefine (component) parameters ================================= */

int nans_match(double a, double b){
  return (*(uint64_t*)&a == *(uint64_t*)&b);
}
int is_unset(double x){
  return nans_match(x, UNSET);
}
int is_set(double x){
  return !nans_match(x, UNSET);
}
int is_valid(double x){
  return !isnan(x)||is_unset(x);
}
int all_unset(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=1;
  for (int i=0; i<n; ++i) if(is_set(va_arg(ptr, double))) ret=0;
  va_end(ptr);
  return ret;
}
int all_set(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=1;
  for (int i=0; i<n; ++i) if(is_unset(va_arg(ptr, double))) ret=0;
  va_end(ptr);
  return ret;
}
int any_unset(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=0;
  for (int i=0; i<n; ++i) if(is_unset(va_arg(ptr, double))) ret=1;
  va_end(ptr);
  return ret;
}
int any_set(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=0;
  for (int i=0; i<n; ++i) if(is_set(va_arg(ptr, double))) ret=1;
  va_end(ptr);
  return ret;
}


/* SECTION: Dynamic Arrays ================================================== */
IArray1d create_iarr1d(int n){
  IArray1d arr2d;
  arr2d = calloc(n, sizeof(int));
  return arr2d;
}
void destroy_iarr1d(IArray1d a){
  free(a);
}

IArray2d create_iarr2d(int nx, int ny){
  IArray2d arr2d;
  arr2d = calloc(nx, sizeof(int *));

  int *p1;
  p1 = calloc(nx*ny, sizeof(int));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}
void destroy_iarr2d(IArray2d a){
  free(a[0]);
  free(a);
}

IArray3d create_iarr3d(int nx, int ny, int nz){
  IArray3d arr3d;
  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(int **));

  // d2
  int **p1;
  p1 = calloc(nx*ny, sizeof(int *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  int *p2;
  p2 = calloc(nx*ny*nz, sizeof(int));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}

void destroy_iarr3d(IArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}

DArray1d create_darr1d(int n){
  DArray1d arr2d;
  arr2d = calloc(n, sizeof(double));
  return arr2d;
}

void destroy_darr1d(DArray1d a){
  free(a);
}

DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}

void destroy_darr2d(DArray2d a){
  free(a[0]);
  free(a);
}

DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;

  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(double **));

  // d2
  double **p1;
  p1 = calloc(nx*ny, sizeof(double *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  double *p2;
  p2 = calloc(nx*ny*nz, sizeof(double));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}

void destroy_darr3d(DArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}


/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Allreduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }, {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

#ifndef STRACPY
/* this is a replacement to strncpy, but ensures that the copy ends with NULL */
/* http://stracpy.blogspot.fr/2011/04/stracpy-strncpy-replacement.html */
#define STRACPY
char *stracpy(char *destination, const char *source, size_t amount)
{
        if (!destination || !source || !amount) return(NULL);
        while(amount--)
          if((*destination++ = *source++) == '\0') break;
        *destination = '\0';
        return destination;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=dirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = dirname ? strlen(dirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, dirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within dirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);

  mem  = mcfull_file(name, ext); /* create dirname/name.ext */

  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;

  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);

    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+");

  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n",
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: detector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m)
    return(detector);

  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (detector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];

  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m)
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin);
      else x = 0;
      if (detector.n && detector.p)
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin);
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw)
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */

      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }

      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {

    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH,
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }

  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n",
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n",
      detector.component, strlen(detector.filename) ? detector.filename : "");

  return(detector);

} /* mcdetector_statistics */

/*******************************************************************************
* detector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=siminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR detector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", instrument_name, instrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=labs(m); n=labs(n); p=labs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, instrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", instrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }


  return(detector);
} /* detector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: siminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < numipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    if (strlen(Parameters) + strlen(ThisParam) + 1 >= CHAR_BUF_LENGTH) break;
    strcat(Parameters, ThisParam);
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, dirname, MC_PATHSEP_C, siminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, instrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);

  fprintf(f, "%sTrace_enabled: %s\n", pre, traceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, defaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre,
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out: output simulation tags/info (both in SIM and data files)
* Used in: siminfo_init (ascii case), mcdetector_out_xD_ascii
*******************************************************************************/
static void mcruninfo_out(char *pre, FILE *f)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre,
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, instrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, dirname ? dirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  // TODO Consider replacing this by a a call to `mcparameterinfo_out(pre+"Param: ", f)`
  /* output parameter string ================================================ */
  for(i = 0; i < numipar; i++) {
      if (mcinputtable[i].par){
	/* Parameters with a default value */
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
	  (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
	  fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        /* ... and those without */
	}else{
	  fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
  fflush(f);
} /* mcruninfo_out */

/*******************************************************************************
 * @brief Print parameter information to the specified file
 * @param pre any beginning-of-line padding
 * @param f the output file
 */
static void mcparameterinfo_out(char * pre, FILE *f){
  if (!f || mcdisable_output_files) return;

  unsigned int nchar = 4;
  for (int i=0; i < numipar; ++i){
    if (mcinputtable[i].par && mcinputtable[i].val && strlen(mcinputtable[i].val) > nchar)
      nchar = strlen(mcinputtable[i].val);
  }
  char * buffer = calloc(nchar+1, sizeof(char));

  if (!buffer) {
    exit(1);
  }

  for (int i=0; i < numipar; ++i) {
    if (mcinputtable[i].par) {
      char * name = mcinputtable[i].name;
      if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
        mcinputtypes[mcinputtable[i].type].printer(buffer, mcinputtable[i].par);
      } else {
        strcpy(buffer, "NULL");
      }
      if (strlen(mcinputtable[i].unit)){
        //fprintf(f, "%s%s %s (\"%s\") = %s\n", pre, mcinputtypes[mcinputtable[i].type].parminfo(name), name, mcinputtable[i].unit, buffer);
        fprintf(f, "%s%s %s/\"%s\" = %s\n", pre, mcinputtypes[mcinputtable[i].type].parminfo(name), name, mcinputtable[i].unit, buffer);
      } else {
        fprintf(f, "%s%s %s = %s\n", pre, mcinputtypes[mcinputtable[i].type].parminfo(name), name, buffer);
      }
    }
  }

  free(buffer);
}

/*******************************************************************************
* siminfo_out:    wrapper to fprintf(siminfo_file)
*******************************************************************************/
void siminfo_out(char *format, ...)
{
  va_list ap;

  if(siminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(siminfo_file, format, ap);
    va_end(ap);
  }
} /* siminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;

  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" :
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f,
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" :
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre,
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);

  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  /* Write data set information to simulation description file. */
  MPI_MASTER(
    siminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", siminfo_file, detector);
    siminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n",
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    siminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", siminfo_file, detector);
    siminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);

      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);

}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        siminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", siminfo_file, detector);
        siminfo_out("end data\n");

        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
        fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      }
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1,
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2,
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0,
          outfile, detector.istransposed);
      }
      fclose(outfile);

      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1,
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=CHAR_BUF_LENGTH; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;

  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);

  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);

    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }

  return(ret);

} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: siminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, instrument_name);

  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK)
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {

    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s",
      dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name,
      instrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   instrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < numipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }

      nxprintattr(f, "name",          instrument_name);
      nxprintf   (f, "name",          instrument_name);
      nxprintattr(f, "Source",        instrument_source);

      nxprintattr(f, "Trace_enabled", traceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  defaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );

      /* add instrument source code when available */
      buffer = mcinfo_readfile(instrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", instrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", instrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)",
          instrument_source, instrument_name);

      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",instrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", instrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name);

      nxprintf   (f, "name",      "%s",     siminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",instrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", dirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif

      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < numipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */

      NXclosegroup(f); /* simulation */
    } /* NXsimulation */

    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[CHAR_BUF_LENGTH];
  if (!f || !detector.m || mcdisable_output_files) return;

  strcpy_valid(data_name,
    detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (siminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */

    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {

      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" :
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" :
                 "xylimits", detector.limits);
      nxprintattr(f, "variables",
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");

      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */

} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[CHAR_BUF_LENGTH];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMPRESSION, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);

    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{

  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int fulldims[3]={detector.m,detector.n,detector.p};
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;

  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);

  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) fulldims[0] = NX_UNLIMITED;

  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, fulldims, NX_COMPRESSION, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */

  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",
        strlen(detector.filename) ? detector.filename : detector.component);
    else
      printf("Append:   \"%s\"\n",
	     strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }

  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);

  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[CHAR_BUF_LENGTH];

  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);

  strcpy_valid(data_name,
    detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (siminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {

      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar,
          1, detector.m, detector.xmin, detector.xmax);

        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar,
          2, detector.n, detector.ymin, detector.ymax);

        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar,
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */

      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);

      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */

  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may.
*   Then the number of recv/send is not constant along nodes, and simulation stalls.
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );

  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];

	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );

  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );

#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* siminfo_init:   open SIM and write header
*******************************************************************************/
FILE *siminfo_init(FILE *f)
{
  int exists=0;

  /* check format */
  if (!mcformat || !strlen(mcformat)
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE")
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }

  /* open the SIM file if not defined yet */
  if (siminfo_file || mcdisable_output_files)
    return (siminfo_file);

#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  siminfo_file = mcnew_file(siminfo_name, "h5", &exists);
    if(!siminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      siminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(siminfo_file); /* points to nxhandle */
  }
#endif

  /* write main description file (only MASTER) */
  MPI_MASTER(

  siminfo_file = mcnew_file(siminfo_name, "sim", &exists);
  if(!siminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    siminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    siminfo_out("%s simulation description file for %s.\n",
      MCCODE_NAME, instrument_name);
    siminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    siminfo_out("Program: %s\n\n", MCCODE_STRING);

    siminfo_out("begin instrument: %s\n", instrument_name);
    mcinfo_out(   "  ", siminfo_file);
    siminfo_out("end instrument\n");

    siminfo_out("\nbegin simulation: %s\n", dirname);
    mcruninfo_out("  ", siminfo_file);
    siminfo_out("end simulation\n");

  }
  return (siminfo_file);

  ); /* MPI_MASTER */

} /* siminfo_init */

/*******************************************************************************
*   siminfo_close:  close SIM
*******************************************************************************/
void siminfo_close()
{
#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif
  if(siminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else {
#endif
      fclose(siminfo_file);
#ifdef USE_NEXUS
    }
#endif
#ifdef USE_MPI
  }
#endif
    siminfo_file = NULL;
  }
} /* siminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = detector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));

} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*
*   t:    title
*   xl:   x-label
*   yl:   y-label
*   xvar: measured variable length
*   x1:   x axus min
*   x2:   x axis max
*   n:    1d data vector lenght
*   p0:   pntr to start of data block#0
*   p1:   pntr to start of data block#1
*   p2:   pntr to start of data block#2
*   f:    filename
*
*   Not included in the macro, and here forwarded to detector_import:
*   c:    ?
*   posa: ?
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  // detector_import calls mcdetector_statistics, which will return different
  // MCDETECTOR versions for 1-D data based on the value of mcformat.
  //
  MCDETECTOR detector = detector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    detector = mcdetector_out_1D_nexus(detector);
  else
#endif
    detector = mcdetector_out_1D_ascii(detector);
  if (detector.p1 != p1 && detector.p1) {
    // mcdetector_statistics allocated memory but it hasn't been freed.
    free(detector.p1);
    // plus undo the other damage done there:
    detector.p0 = p0; // was set to NULL
    detector.p1 = p1; // was set to this_p1
    detector.p2 = p2; // was set to NULL
    detector.m = detector.n; // (e.g., labs(n))
    detector.n = 1;  // not (n x n)
    detector.istransposed = n < 0 ? 1 : 0;
  }
  return detector;

} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   Special case for list: master creates file first, then slaves append their
*   blocks without header-
*
*   t:    title
*   xl:   x-label
*   yl:   y-label
*   x1:   x axus min
*   x2:   x axis max
*   y1:   y axis min
*   y2:   y axis max
*   m:    dim 1 (x) size
*   n:    dim 2 (y) size
*   p0:   pntr to start of data block#0
*   p1:   pntr to start of data block#1
*   p2:   pntr to start of data block#2
*   f:    filename
*
*   Not included in the macro, and here forwarded to detector_import:
*   c:    ?
*   posa: ?
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];

  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (labs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (labs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));

} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;

  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,labs(m),1,labs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);

  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    dirname = dir;
  else
    dirname = dir+strlen("file://");


#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif
    if(mkdir(dirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
#ifdef USE_MPI
    }
#endif

  /* remove trailing PATHSEP (if any) */
  while (strlen(dirname) && dirname[strlen(dirname) - 1] == MC_PATHSEP_C)
    dirname[strlen(dirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", instrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", dirname ? dirname : ".");
  mcruninfo_out("  ", stdout);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

/*******************************************************************************
* mcparameterinfo: display instrument parameter info to stdout and exit
*******************************************************************************/
static void
mcparameterinfo(void)
{
  mcparameterinfo_out("  ", stdout);
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcparameterinfo */



#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays */
/* Within the TRACE scope we are now using _particle->uid directly */
unsigned long long int mcget_run_num() // shuld be (_class_particle* _particle) somehow
{
  /* This function only remains for the few cases outside TRACE where we need to know
     the number of simulated particles */
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(!mcseed) {
  //  srandom(mcseed);
  //} else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  // Do nothing here, better use interactive zoom from the tools
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* Draws a circle with center (x,y,z), radius (r), and in the plane
 * with normal (nx,ny,nz)*/
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz){
    int i;
    if(nx==0 && ny && nz==0){
        for (i=0;i<24; i++){
            mcdis_line(x+r*sin(i*2*PI/24),y,z+r*cos(i*2*PI/24),
                    x+r*sin((i+1)*2*PI/24),y,z+r*cos((i+1)*2*PI/24));
        }
    }else{
        double mx,my,mz;
        /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
        /*draw circle*/
        for (i=0;i<24; i++){
            double ux,uy,uz;
            double wx,wy,wz;
            rotate(ux,uy,uz, mx,my,mz, i*2*PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*PI/24, nx,ny,nz);
            mcdis_line(x+ux*r,y+uy*r,z+uz*r,
                    x+wx*r,y+wy*r,z+wz*r);
        }
    }
}

/* Draws a cylinder with center at (x,y,z) with extent (r,height).
 * The cylinder axis is along the vector nx,ny,nz.
 * N determines how many vertical lines are drawn.*/
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz){
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    NORM(nx,ny,nz);
    double h_2=height/2.0;
    mcdis_Circle(x+nx*h_2,y+ny*h_2,z+nz*h_2,r,nx,ny,nz);
    mcdis_Circle(x-nx*h_2,y-ny*h_2,z-nz*h_2,r,nx,ny,nz);

    double mx,my,mz;
    /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
    if(nx==0 && ny && nz==0){
        mx=my=0;mz=1;
    }else{
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
    }
    /*draw circle*/
    for (i=0; i<24; i++){
        double ux,uy,uz;
        rotate(ux,uy,uz, mx,my,mz, i*2*PI/24, nx,ny,nz);
        mcdis_line(x+nx*h_2+ux*r, y+ny*h_2+uy*r, z+nz*h_2+uz*r,
                 x-nx*h_2+ux*r, y-ny*h_2+uy*r, z-nz*h_2+uz*r);
    }
}

/* draws a sphere with center at (x,y,z) with extent (r)
 * The sphere is drawn using N longitudes and N latitudes.*/
void mcdis_sphere(double x, double y, double z, double r, int N){
    double nx,ny,nz;
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    nx=0;ny=0;nz=1;
    mcdis_Circle(x,y,z,r,nx,ny,nz);
    for (i=1;i<N;i++){
        rotate(nx,ny,nz, nx,ny,nz, PI/N, 0,1,0);
        mcdis_Circle(x,y,z,r,nx,ny,nz);
    }
    /*lastly draw a great circle perpendicular to all N circles*/
    //mcdis_Circle(x,y,z,radius,1,0,0);

    for (i=1;i<=N;i++){
        double yy=-r+ 2*r*((double)i/(N+1));
        mcdis_Circle(x,y+yy ,z,  sqrt(r*r-yy*yy) ,0,1,0);
    }
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {
  #ifndef OPENACC
  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  #endif
  return;
}

mcstatic void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/* coords_test_zero: check if zero vector*/
int coords_test_zero(Coords a){
  return ( a.x==0 && a.y==0 && a.z==0 );
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}


/* SECTION: GPU algorithms ================================================== */


/*
*  Divide-and-conquer strategy for parallelizing this task: Sort absorbed
*  particles last.
*
*   particles:  the particle array, required to checking _absorbed
*   pbuffer:    same-size particle buffer array required for parallel sort
*   len:        sorting area-of-interest size (e.g. from previous calls)
*   buffer_len: total array size
*   flag_split: if set, multiply live particles into absorbed slots, up to buffer_len
*   multiplier: output arg, becomes the  SPLIT multiplier if flag_split is set
*/
#ifdef FUNNEL
long sort_absorb_last(_class_particle* particles, _class_particle* pbuffer, long len, long buffer_len, long flag_split, long* multiplier) {
  #define SAL_THREADS 1024 // num parallel sections
  if (len<SAL_THREADS) return sort_absorb_last_serial(particles, len);

  if (multiplier != NULL) *multiplier = -1; // set default out value for multiplier
  long newlen = 0;
  long los[SAL_THREADS]; // target array startidxs
  long lens[SAL_THREADS]; // target array sublens
  long l = floor(len/(SAL_THREADS-1)); // subproblem_len
  long ll = len - l*(SAL_THREADS-1); // last_subproblem_len

  // TODO: The l vs ll is too simplistic, since ll can become much larger
  // than l, resulting in idling. We should distribute lengths more evenly.

  // step 1: sort sub-arrays
  #pragma acc parallel loop present(particles, pbuffer)
  for (unsigned long tidx=0; tidx<SAL_THREADS; tidx++) {
    long lo = l*tidx;
    long loclen = l;
    if (tidx==(SAL_THREADS-1)) loclen = ll; // last sub-problem special case
    long i = lo;
    long j = lo + loclen - 1;

    // write into pbuffer at i and j
    #pragma acc loop seq
    while (i < j) {
      #pragma acc loop seq
      while (!particles[i]._absorbed && i<j) {
        pbuffer[i] = particles[i];
        i++;
      }
      #pragma acc loop seq
      while (particles[j]._absorbed && i<j) {
        pbuffer[j] = particles[j];
        j--;
      }
      if (i < j) {
        pbuffer[j] = particles[i];
        pbuffer[i] = particles[j];
        i++;
        j--;
      }
    }
    // transfer edge case
    if (i==j)
      pbuffer[i] = particles[i];

    lens[tidx] = i - lo;
    if (i==j && !particles[i]._absorbed) lens[tidx]++;
  }

  // determine lo's
  long accumlen = 0;
  #pragma acc loop seq
  for (long idx=0; idx<SAL_THREADS; idx++) {
    los[idx] = accumlen;
    accumlen = accumlen + lens[idx];
  }

  // step 2: write non-absorbed sub-arrays to psorted/output from the left
  #pragma acc parallel loop present(pbuffer)
  for (unsigned long tidx=0; tidx<SAL_THREADS; tidx++) {
    long j, k;
    #pragma acc loop seq
    for (long i=0; i<lens[tidx]; i++) {
      j = i + l*tidx;
      k = i + los[tidx];
      particles[k] = pbuffer[j];
    }
  }
  //for (int ii=0;ii<accumlen;ii++) printf("%ld ", (psorted[ii]->_absorbed));

  // return (no SPLIT)
  if (flag_split != 1)
    return accumlen;

  // SPLIT - repeat the non-absorbed block N-1 times, where len % accumlen = N + R
  int mult = buffer_len / accumlen; // TODO: possibly use a new arg, bufferlen, rather than len

  // not enough space for full-block split, return
  if (mult <= 1)
    return accumlen;

  // copy non-absorbed block
  #pragma acc parallel loop present(particles)
  for (long tidx = 0; tidx < accumlen; tidx++) { // tidx: thread index
    unsigned long randstate[7];
    _class_particle sourcebuffer;
    _class_particle targetbuffer;
    // assign reduced weight to all particles
    particles[tidx].p=particles[tidx].p/mult;
    #pragma acc loop seq
    for (long bidx = 1; bidx < mult; bidx++) { // bidx: block index
      // preserve absorbed particle (for randstate)
      sourcebuffer = particles[bidx*accumlen + tidx];
      // buffer full particle struct
      targetbuffer = particles[tidx];
      // reassign previous randstate
      targetbuffer.randstate[0] = sourcebuffer.randstate[0];
      targetbuffer.randstate[1] = sourcebuffer.randstate[1];
      targetbuffer.randstate[2] = sourcebuffer.randstate[2];
      targetbuffer.randstate[3] = sourcebuffer.randstate[3];
      targetbuffer.randstate[4] = sourcebuffer.randstate[4];
      targetbuffer.randstate[5] = sourcebuffer.randstate[5];
      targetbuffer.randstate[6] = sourcebuffer.randstate[6];
      // apply
      particles[bidx*accumlen + tidx] = targetbuffer;
    }
  }

  // set out split multiplier value
  *multiplier = mult;

  // return expanded array size
  return accumlen * mult;
}

#endif

/*
*  Fallback serial version of the one above.
*/
long sort_absorb_last_serial(_class_particle* particles, long len) {
  long i = 0;
  long j = len - 1;
  _class_particle pbuffer;

  // bubble
  while (i < j) {
    while (!particles[i]._absorbed && i<j) i++;
    while (particles[j]._absorbed && i<j) j--;
    if (i < j) {
      pbuffer = particles[j];
      particles[j] = particles[i];
      particles[i] = pbuffer;
      i++;
      j--;
    }
  }

  // return new length
  if (i==j && !particles[i]._absorbed)
    return i + 1;
  else
    return i;
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void mccoordschange(Coords a, Rotation t, _class_particle *particle)
{
  Coords b, c;

  b.x = particle->x;
  b.y = particle->y;
  b.z = particle->z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  particle->x = b.x;
  particle->y = b.y;
  particle->z = b.z;

#if MCCODE_PARTICLE_CODE == 2112
    if (particle->vz != 0.0 || particle->vx != 0.0 || particle->vy != 0.0)
      mccoordschange_polarisation(t, &(particle->vx), &(particle->vy), &(particle->vz));

    if (particle->sz != 0.0 || particle->sx != 0.0 || particle->sy != 0.0)
      mccoordschange_polarisation(t, &(particle->sx), &(particle->sy), &(particle->sz));
#elif MCCODE_PARTICLE_CODE == 22
    if (particle->kz != 0.0 || particle->kx != 0.0 || particle->ky != 0.0)
      mccoordschange_polarisation(t, &(particle->kx), &(particle->ky), &(particle->kz));

    if (particle->Ez != 0.0 || particle->Ex != 0.0 || particle->Ey != 0.0)
      mccoordschange_polarisation(t, &(particle->Ex), &(particle->Ey), &(particle->Ez));
#endif
}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
void normal_vec(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order_old(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

int solve_2nd_order(double *t0, double *t1, double A, double B, double C){
  int retval=0;
  double sign=copysign(1.0,B);
  double dt0,dt1;

  dt0=0;
  dt1=0;
  if(t1){ *t1=0;}

  /*protect against rounding errors by locally equating DBL_EPSILON with 0*/
  if (fabs(A)<DBL_EPSILON){
    A=0;
  }
  if (fabs(B)<DBL_EPSILON){
    B=0;
  }
  if (fabs(C)<DBL_EPSILON){
    C=0;
  }

  /*check if coefficient are sane*/
  if( A==0  && B==0){
    retval=0;
  }else{
    if(A==0){
      /*equation is linear*/
      dt0=-C/B;
      retval=1;
    }else if (C==0){
      /*one root is 0*/
      if(sign<0){
        dt0=0;dt1=-B/A;
      }else{
        dt0=-B/A;dt1=0;
      }
      retval=2;
    }else{
      /*a regular 2nd order eq. Also works out fine for B==0.*/
      double D;
      D=B*B-4*A*C;
      if (D>=0){
        dt0=(-B - sign*sqrt(B*B-4*A*C))/(2*A);
        dt1=C/(A*dt0);
        retval=2;
      }else{
        /*no real roots*/
        retval=0;
      }
    }
    /*sort the solutions*/
    if (retval==1){
      /*put both solutions in t0 and t1*/
      *t0=dt0;
      if(t1) *t1=dt1;
    }else{
      /*we have two solutions*/
      /*swap if both are positive and t1 smaller than t0 or t1 the only positive*/
      int swap=0;
      if(dt1>0 && ( dt1<dt0 || dt0<=0) ){
        swap=1;
      }
      if (swap){
        *t0=dt1;
        if(t1) *t1=dt0;
      }else{
        *t0=dt0;
        if(t1) *t1=dt0;
      }
    }

  }
  return retval;

} /*solve_2nd_order_improved*/


/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void _randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double radius,
        _class_particle* _particle)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos(1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
}
/* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void _randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double width, double height, Rotation A,
        _class_particle* _particle)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);
}
/* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/
void _randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi,
        double width, double height, Rotation A,
        double lx, double ly, double lz, int order,
        _class_particle* _particle)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
        *solid_angle = *solid_angle * cos_theta;
      }
    }
  }
}
/* randvec_target_rect_real */


/* SECTION: random numbers ==================================================

  How to add a new RNG:

  - Use an rng with a manegable state vector, e.g. of lengt 4 or 7. The state
  will sit on the particle struct as a "randstate_t state[RANDSTATE_LEN]"
  - If the rng has a long state (as MT), set an empty "srandom" and initialize
  it explicitly using the appropriate define (RNG_ALG)
  - Add a seed and a random function (the transforms will be reused)
  - Write the proper defines in mccode-r.h, e.g. randstate_t and RANDSTATE_LEN,
  srandom and random.
  - Compile using -DRNG_ALG=<selector int value>

============================================================================= */


/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */
/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/
#include <stdio.h>
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

unsigned long mt[N]; /* the array for the state vector  */
int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

// required for common rng alg interface (see RNG_ALG usage in mccode-r.h)
void mt_srandom_empty() {}

// initializes mt[N] with a seed
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}
/* Initialize by an array with array-length.
   Init_key is the array for initializing keys.
   key_length is its length. */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}
/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK
/* End of "Mersenne Twister". */


/*
KISS

 From: http://www.helsbreth.org/random/rng_kiss.html
 Scott Nelson 1999

 Based on Marsaglia's KISS or (KISS+SWB) <http://www.cs.yorku.ca/~oz/marsaglia-
rng.html>

 KISS - Keep it Simple Stupid PRNG

 the idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
        x(n)=a*x(n-1)+1 mod 2^32
        y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
        z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
      2^32*(2^32-1)*(2^63+2^32-1) > 2^127
*/

/* the KISS state is stored as a vector of 7 unsigned long        */
/*   0  1  2  3  4      5  6   */
/* [ x, y, z, w, carry, k, m ] */

unsigned long *kiss_srandom(unsigned long state[7], unsigned long seed) {
  if (seed == 0) seed = 1;
  state[0] = seed | 1; // x
  state[1] = seed | 2; // y
  state[2] = seed | 4; // z
  state[3] = seed | 8; // w
  state[4] = 0;        // carry
  return 0;
}

unsigned long kiss_random(unsigned long state[7]) {
    state[0] = state[0] * 69069 + 1;
    state[1] ^= state[1] << 13;
    state[1] ^= state[1] >> 17;
    state[1] ^= state[1] << 5;
    state[5] = (state[2] >> 2) + (state[3] >> 3) + (state[4] >> 2);
    state[6] = state[3] + state[3] + state[2] + state[4];
    state[2] = state[3];
    state[3] = state[6];
    state[4] = state[5] >> 30;
    return state[0] + state[1] + state[3];
}
/* end of "KISS" rng */


/* FAST KISS in another implementation (Hundt) */

//////////////////////////////////////////////////////////////////////////////
// fast keep it simple stupid generator
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Thomas Mueller hash for initialization of rngs
// http://stackoverflow.com/questions/664014/
//        what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
//////////////////////////////////////////////////////////////////////////////
randstate_t _hash(randstate_t x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x);
  return x;
}


// SECTION: random number transforms ==========================================



// generate a random number from normal law
double _randnorm(randstate_t* state)
{
  static double v1, v2, s; /* removing static breaks comparison with McStas <= 2.5 */
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = _rand01(state);
      u2 = _rand01(state);
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}
// another one
double _randnorm2(randstate_t* state) {
  double x, y, r;
  do {
      x = 2.0 * _rand01(state) - 1.0;
      y = 2.0 * _rand01(state) - 1.0;
      r = x*x + y*y;
  } while (r == 0.0 || r >= 1.0);
  return x * sqrt((-2.0 * log(r)) / r);
}

// Generate a random number from -1 to 1 with triangle distribution
double _randtriangle(randstate_t* state) {
	double randnum = _rand01(state);
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}
double _rand01(randstate_t* state) {
	double randnum;
	randnum = (double) _random();
  // TODO: can we mult instead of div?
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}
// Return a random number between 1 and -1
double _randpm1(randstate_t* state) {
	double randnum;
	randnum = (double) _random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}
// Return a random number between 0 and max.
double _rand0max(double max, randstate_t* state) {
	double randnum;
	randnum = (double) _random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}
// Return a random number between min and max.
double _randminmax(double min, double max, randstate_t* state) {
	return _rand0max(max - min, state) + max;
}


/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", instrument_name, instrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of particles to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --list-parameters          Print the instrument parameters to standard out\n"
"  --meta-list                Print names of components which defined metadata\n"
"  --meta-defined COMP[:NAME] Print component defined metadata names, or (0,1) if NAME provided\n"
"  --meta-type COMP:NAME      Print metadata format type specified in definition\n"
"  --meta-data COMP:NAME      Print the metadata text\n"
"  --source                   Show the instrument code which was compiled.\n"
#ifdef OPENACC
"\n"
"  --vecsize                  OpenACC vector-size (default: 128)\n"
"  --numgangs                 Number of OpenACC gangs (default: 7813)\n"
"  --gpu_innerloop            Maximum rays to process pr. OpenACC \n"
"                             kernel run (default: 2147483647)\n"
"\n"
#endif
"\n"
"  --bufsiz                   Monitor_nD list/buffer-size (default: 1000000)\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
#ifdef OPENACC
  fprintf(stderr,
  "This instrument has been compiled with NVIDIA GPU support through OpenACC.\n  Running on systems without such devices will lead to segfaults.\nFurter, fprintf, sprintf and printf have been removed from any component TRACE.\n");
#endif

  if(numipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < numipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(traceenabled) {
  mcdotrace = 1;
  #pragma acc update device ( mcdotrace )
 } else {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    instrument_name, instrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to numipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((numipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < numipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]) || !strcmp("--version", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if (!strcmp("--list-parameters", argv[i]))
      mcparameterinfo();
    else if (!strcmp("--meta-list", argv[i]) && ((i+1) >= argc || argv[i+1][0] == '-')){
      //printf("Components with metadata defined:\n");
      exit(metadata_table_print_all_components(num_metadata, metadata_table) == 0);
    }
    else if (!strcmp("--meta-defined", argv[i]) && (i+1) < argc){
      exit(metadata_table_print_component_keys(num_metadata, metadata_table, argv[i+1]) == 0);
    }
    else if (!strcmp("--meta-type", argv[i]) && (i+1) < argc){
      char * literal_type = metadata_table_type(num_metadata, metadata_table, argv[i+1]);
      if (literal_type == NULL) exit(1);
      printf("%s\n", literal_type);
      exit(0);
    }
    else if (!strcmp("--meta-data", argv[i]) && (i+1) < argc){
      char * literal = metadata_table_literal(num_metadata, metadata_table, argv[i+1]);
      if (literal == NULL) exit(1);
      printf("%s\n", literal);
      exit(0);
    }
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]) || !strcmp("--verbose", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strncmp("--vecsize=", argv[i], 10)) {
      vecsize=atoi(&argv[i][10]);
    }    
    else if(!strcmp("--vecsize", argv[i]) && (i + 1) < argc) {
      vecsize=atoi(argv[++i]);
    }
    else if(!strcmp("--bufsiz", argv[i]) && (i + 1) < argc) {
      MONND_BUFSIZ=atoi(argv[++i]);
    }
    else if(!strncmp("--numgangs=", argv[i], 11)) {
      numgangs=atoi(&argv[i][11]);
    }
    else if(!strcmp("--numgangs", argv[i]) && (i + 1) < argc) {
      numgangs=atoi(argv[++i]);
    }
    else if(!strncmp("--gpu_innerloop=", argv[i], 16)) {
      gpu_innerloop=(long)strtod(&argv[i][16], NULL);
    }
    else if(!strcmp("--gpu_innerloop", argv[i]) && (i + 1) < argc) {
      gpu_innerloop=(long)strtod(argv[++i], NULL);
    }

    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(!strcmp("--source", argv[i])) {
      printf("/* Source code %s from %s: */\n"
        "/******************************************************************************/\n"
        "%s\n"
        "/******************************************************************************/\n"
        "/* End of source code %s from %s */\n",
        instrument_name, instrument_source, instrument_code,
        instrument_name, instrument_source);
      exit(1);
    }
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < numipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == numipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < numipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir) && !mcdisable_output_files) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", instrument_name, instrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    save(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    finally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  init();

  /* *** parse options *** */
  SIG_MESSAGE("[" __FILE__ "] main START");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

/* embedding file "mcstas-r.c" */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables
*******************************************************************************/
_class_particle mcsetstate(double x, double y, double z, double vx, double vy, double vz,
			   double t, double sx, double sy, double sz, double p, int mcgravitation, void *mcMagnet, int mcallowbackprop)
{
  _class_particle mcneutron;

  mcneutron.x  = x;
  mcneutron.y  = y;
  mcneutron.z  = z;
  mcneutron.vx = vx;
  mcneutron.vy = vy;
  mcneutron.vz = vz;
  mcneutron.t  = t;
  mcneutron.sx = sx;
  mcneutron.sy = sy;
  mcneutron.sz = sz;
  mcneutron.p  = p;
  mcneutron.mcgravitation = mcgravitation;
  mcneutron.mcMagnet = mcMagnet;
  mcneutron.allow_backprop = mcallowbackprop;
  mcneutron._uid       = 0;
  mcneutron._index     = 1;
  mcneutron._absorbed  = 0;
  mcneutron._restore   = 0;
  mcneutron._scattered = 0;

  return(mcneutron);
} /* mcsetstate */

/*******************************************************************************
* mcgetstate: get neutron parameters from particle structure
*******************************************************************************/
_class_particle mcgetstate(_class_particle mcneutron, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
  *x  =  mcneutron.x;
  *y  =  mcneutron.y;
  *z  =  mcneutron.z;
  *vx =  mcneutron.vx;
  *vy =  mcneutron.vy;
  *vz =  mcneutron.vz;
  *t  =  mcneutron.t;
  *sx =  mcneutron.sx;
  *sy =  mcneutron.sy;
  *sz =  mcneutron.sz;
  *p  =  mcneutron.p;

  return(mcneutron);
} /* mcgetstate */


/*******************************************************************************
* mcgenstate: set default neutron parameters
*******************************************************************************/
// Moved to generated code
/* #pragma acc routine seq */
/* _class_particle mcgenstate(void) */
/* { */
/*   return(mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, mcgravitation, mcMagnet, mcallowbackprop)); */
/* } */

/*******************************************************************************
* mccoordschanges: old style rotation routine rot -> (x y z) ,(vx vy vz),(sx,sy,sz)
*******************************************************************************/
void
mccoordschanges(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) )
    mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) )
    mccoordschange_polarisation(t, sx, sy, sz);

}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight)
* return 0 if outside and 1 if inside
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999.
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98
  *******************************************************************************/
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1
 *******************************************************************************/
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */


/* *****************************************************************************
* Start of instrument 'template_simple' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "C:\\mcstas-3.4\\lib\\"
int   defaultmain         = 1;
char  instrument_name[]   = "template_simple";
char  instrument_source[] = "WARP.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument template_simple source code WARP.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'template_simple' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM focus;
  MCNUM E_m;
  MCNUM d_E;
  MCNUM thickness;
  MCNUM sample_x;
  MCNUM sample_y;
  MCNUM L0;
  MCNUM L1;
  MCNUM Ld;
  MCNUM det_rot;
  MCNUM slit_x;
  MCNUM slit_y;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT and GROUP control logic */
struct instrument_logic_struct {
  long Group_Analyzer; /* equals index of scattering comp when in group */
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'template_simple' */
/* Counters per component instance */
  double counter_AbsorbProp[57]; /* absorbed events in PROP routines */
  double counter_N[57], counter_P[57], counter_P2[57]; /* event counters after each component instance */
  _class_particle _trajectory[57]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[57]; /* positions of all components */
  Coords _position_absolute[57];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 12;
struct mcinputtable_struct mcinputtable[] = {
  "focus", &(_instrument_var._parameters.focus), instr_type_double, "4", "",
  "E_m", &(_instrument_var._parameters.E_m), instr_type_double, "5.325", "",
  "d_E", &(_instrument_var._parameters.d_E), instr_type_double, "1.375", "",
  "thickness", &(_instrument_var._parameters.thickness), instr_type_double, "0.002", "",
  "sample_x", &(_instrument_var._parameters.sample_x), instr_type_double, "0.005", "",
  "sample_y", &(_instrument_var._parameters.sample_y), instr_type_double, "0.015", "",
  "L0", &(_instrument_var._parameters.L0), instr_type_double, "1.1749", "",
  "L1", &(_instrument_var._parameters.L1), instr_type_double, "2.1988499999999997", "",
  "Ld", &(_instrument_var._parameters.Ld), instr_type_double, "0.5375925315701473", "",
  "det_rot", &(_instrument_var._parameters.det_rot), instr_type_double, "3.0291893730987205", "",
  "slit_x", &(_instrument_var._parameters.slit_x), instr_type_double, "0.04", "",
  "slit_y", &(_instrument_var._parameters.slit_y), instr_type_double, "0.04", "",
  NULL, NULL, instr_type_double, ""
};

struct metadata_table_struct metadata_table[] = {
  "", "", "", ""
};
int num_metadata = 0;

/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'Monochromator_bent'. */
	#include <string.h>

////////////////////////////
// Mathematical Functions //
////////////////////////////

	// Function to return sign of a double
	double sign(double x){
		if (x > 0) return 1;
		if (x < 0) return -1;
		if (x == 0.0) return 0;
	}

	// Function to return a value squared
	double square(double x){
		return x*x;
	}

	// Function to return a random point on a gaussian 
	// with variance sigma and average at 0
	double gauss_0(double sigma){
		double u1, u2;
		u1 = rand01();
		u2 = rand01();
		double r = sqrt(-2 * log(u1));
		double theta = 2 * M_PI * u2;
		return sigma * r * cos(theta);
	}

	// The following two function returns, respectively, the Gaussian cumulative distribution function,
	// And the inverse gaussian cumulative distribution function. The latter is a numeric approximation
	// that could be optimized, but not neccisarily any need for it
	double normalCDF(double value) {
		return 0.5 * erfc(-value * M_SQRT1_2);
	}

	double inverse_of_normal_cdf(double p){
		if (p <= 0 || p >= 1) return sign(p)*6;

		double mu = 0;
		double sigma = 1;
		double r, val;
		double q = p - 0.5;

		if (abs(q) <= .425) {
			r = .180625 - q * q;
			val =
				q * (((((((r * 2509.0809287301226727 +
					33430.575583588128105) * r + 67265.770927008700853) * r +
					45921.953931549871457) * r + 13731.693765509461125) * r +
					1971.5909503065514427) * r + 133.14166789178437745) * r +
					3.387132872796366608)
				/ (((((((r * 5226.495278852854561 +
					28729.085735721942674) * r + 39307.89580009271061) * r +
					21213.794301586595867) * r + 5394.1960214247511077) * r +
					687.1870074920579083) * r + 42.313330701600911252) * r + 1);
		}
		else {
			if (q > 0) {
				r = 1 - p;
			}
			else {
				r = p;
			}

			r = sqrt(-log(r));

			if (r <= 5) 
			{
				r += -1.6;
				val = (((((((r * 7.7454501427834140764e-4 +
					.0227238449892691845833) * r + .24178072517745061177) *
					r + 1.27045825245236838258) * r +
					3.64784832476320460504) * r + 5.7694972214606914055) *
					r + 4.6303378461565452959) * r +
					1.42343711074968357734)
					/ (((((((r *
						1.05075007164441684324e-9 + 5.475938084995344946e-4) *
						r + .0151986665636164571966) * r +
						.14810397642748007459) * r + .68976733498510000455) *
						r + 1.6763848301838038494) * r +
						2.05319162663775882187) * r + 1);
			}
			else { /* very close to  0 or 1 */
				r += -5;
				val = (((((((r * 2.01033439929228813265e-7 +
					2.71155556874348757815e-5) * r +
					.0012426609473880784386) * r + .026532189526576123093) *
					r + .29656057182850489123) * r +
					1.7848265399172913358) * r + 5.4637849111641143699) *
					r + 6.6579046435011037772)
					/ (((((((r *
						2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
						r + 1.8463183175100546818e-5) * r +
						7.868691311456132591e-4) * r + .0148753612908506148525)
						* r + .13692988092273580531) * r +
						.59983220655588793769) * r + 1);
			}

			if (q < 0.0) {
				val = -val;
			}
		}

		return mu + sigma * val;
	}
	// Integral needed for debye factor #TODO: Deeper explaination is needed here.
	double calculate_phi_integral(double x){
		// Asymptotic approximation
		if (x > 5) return PI * PI / 6 - exp(-x)/(x+1);
		double z = 1 + x/(exp(x)-1);
		double dx = x/100;
		double ksi;
		for (int i = 2; i <= 100; i++) {
			ksi = (i-1)*dx;
			switch (i%2){
				case 1:
					z = z + 4 * ksi/(exp(ksi)-1);
					break;
				case 0:
					z = z + 2 * ksi/(exp(ksi)-1);
					break;
			}
		}
		return z*dx/3;
	}

//////////////////////////////////////////////
// Structures for neutron and monochromator //
//////////////////////////////////////////////

	// Neutron values to save throughout the simulation
	struct neutron_values {
								// **Values from McStas, and calculated values related to this**								
		double r[3]; 			// Position at a given time during the simulation
		double v[3]; 			// velocity of neutron
		double v_size; 			// speed
		double l; 				// De Broglie wavelength of neutron
		double E; 				// energy of neutron
		double ki[3]; 			// Incoming wavevector
		double ki_size; 		// size of incoming wavevector

								// **Final wavevector, used to update the final speed**
		double kf[3]; 			// outgoing wavevector
		double kf_size; 		// size of outgoing wavevector
		double Bragg_angle;		// Bragg angle of scattering
		
								// **Calcualted scattering vector, see Saroun eq. 2 & 4**
		double Gr[3]; 			// Total scattering vector 
		double Gr_size;			// Size of scattering vector
		double Grad_G[3];		// Gradient component of scattering vector
		double Grad_G_size;		// Size of gradient component
		double perp_to_Gk[3];	// Vector perpendicular to G and ki
		double perp_size;		// Normalizing factor for perpendicular vector
		double gamma_0;			// Mosaicity in perpendicular direction
		double Ggamma[3];		// Perpendicular mosaicity
		double eps;				// Mosaicity in plane with G and ki
		double eps_zero;		// 0th order expansion of epsilon
		double beta;			// 1st order expansion of epsilon with respect to path travelled through crystal

								// **Path, time and probability of neutron through monochromator, see Saroun eq. 5 & 8**
		double path; 			// Length of the path the neutron follows
		double max_tau;			// Maximal possible time through crystal path
		double entry_time; 		// Time to start of crystal
		double exit_time; 		// Time to end of crystal (if not scattered)
		double tau[5];			// Time traveled each reflection. 
		double P[5];			// Probability of each reflection
		int scatter_count;	// Number of reflections

		int direction;			// Direction travelling (1 if towards monochromator, -1 if away)
	};		

	// Monochromator Values to save throughout the simulation
	struct monochromator_values{
													// **Set by user**
		double length, height, thickness;			// Dimensions of Monochromator
		double mos;									// Mosaicity
		int type;									// Value to define if flat, bent, and/or mosaic
		double radius_horizontal;					// Radius of curvature
		double temperature_mono;					// Temperature of monochromator
		double domain_thickness;					// Thickness of crystal domains

													// **Defined by Material, from RESTRAX**
		double Debye_Waller_factor;					// Debye Waller Factor for reflectivity
		double lattice_spacing;						// Lattice Spacing
		double Maier_Leibnitz_reflectivity;			// Maier Leibnitz Reflectivity							
		double poisson_ratio;						// Poisson ratio for curvature of material
		double bound_atom_scattering_cross_section;	// Scattering cross section - bound						
		double absorption_for_1AA_Neutrons;			// Scattering cross section - absportion for 1Å						
		double incoherent_scattering_cross_section;	// Scattering cross section - incoherent								
		double volume;								// Unit cell volume	
		double Constant_from_Freund_paper;			// Constant from Freund paper for different materials as monochroamtors						
		double debye_temperature;					// Debye Temperature				
		double atomic_number;						// Atomic number			
		double B0;									// Debye-Waller factor at 0 K
		double BT;									// Debye-Waller factor at set temperature
		double single_phonon_absorption;			// Scattering cross section - absportion for single phonon					
		double multiple_phonon_absorption;			// Scattering cross section - absportion for multiple phonon											
		double nuclear_capture_absorption;			// Scattering cross section - absportion for nuclear capture											

													// **Scattering vector cosntants**
		double G_size_zero;							// Scattering vector for flat monochromator
		double G0[3];								// Scattering vector tilted by off-cut of material
		double perp_G[3];							// Perpendicular to G0 and ki
		double lattice_spacing_gradient_field[3][3];// Gradient for scattering vector (curvature dependent)						
		double max_angle;							// Check if inside monochromator
		double min_angle;							// Check if inside monochromator	
	};

///////////////////////////////
// Crystal Plane information //
///////////////////////////////

	// All information about which crystal is being simulated.
	// Most of this information is used to calculate the attenuation and curvature effects.
	// Additional crystals cann be added if one wants, but information must be available for 
	// these materials.

	// Enumerated planes
	enum crystal_plane {Cu111, Cu200, Cu220, Cu311, Cu400, Cu331, Cu420, Cu440, Ge111, Ge220, Ge311,
				Ge400, Ge331, Ge422, Ge511, Ge533, Ge711, Ge551, Si111, Si220, Si311, Si400, Si331, 
				Si422, Si333, Si511, Si440, Si711, Si551, Be10, Be100, Be102, Be103, Be110, Be112, Be200, 
				Be00_2, Be10_1, PG00_2,PG00_4,PG00_6, Fe110, HS111,HS222,HS111star,Di111,Di220, Di311, Di400, 
				Di331, Di422, Di333, Di511, Di440};

	// An array containing all the possible strings that will be accepted if given as an 
	// argument to the parameter plane_of_reflection 
	const char* crystal_planeStrings[] = {
		"Cu111", "Cu200", "Cu220", "Cu311", "Cu400", "Cu331", "Cu420", "Cu440", "Ge111", "Ge220", "Ge311",
				"Ge400", "Ge331", "Ge422", "Ge511", "Ge533", "Ge711", "Ge551", "Si111", "Si220", "Si311", "Si400", "Si331", 
				"Si422", "Si333", "Si511", "Si440", "Si711", "Si551"," Be10", "Be100", "Be102", "Be103", "Be110", "Be112", "Be200", 
				"Be00_2", "Be10_1", "PG00_2","PG00_4","PG00_6", "Fe110", "HS111","HS222","HS111star","Di111","Di220", "Di311", "Di400", 
				"Di331", "Di422", "Di333", "Di511", "Di440"};

	// Function to convert a string to an enum value
	enum crystal_plane stringToEnum(const char* plane) {
		for (int i = 0; i < sizeof(crystal_planeStrings) / sizeof(crystal_planeStrings[0]); ++i) {
			if (strcmp(plane, crystal_planeStrings[i]) == 0) {
				return (enum crystal_plane)i;
			}
		}
	}

	/* Crystal table for perfect crystal bent monochromator
		Table copied from SIMRES, current url: https://github.com/saroun/simres
		Contents: dhkl, QML,sigmab,sigmaa,V0,A,thetaD,C2,poi
		dhkl ... Lattice spacing of crystal plane.
		QML = 4*PI*(F*dhkl/V0)**2 [ A^-1 cm^-1]
		sigmab ... bound-atom scattering cross-section [barn]
		sigmaa ... absorption for 1A neutrons [barn*A^-1]
		sigmai ... incoherent scattering cross-section [barn]
		V0 .... volume [A^3]/atom
		A  .... atomic number
		thetaD .... Debye temperature (K)
		C2 .... constant from the Freund's paper  [A^-2 eV^-1]
		poi .... Poisson elastic constant */
	double crystal_table[56][10] = {{ 2.087063,  0.23391E+00 ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 1.80745 , 0.17544E+00  ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 1.27806 , 0.87718E-01  ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 1.089933,  0.63795E-01 ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 0.903725,  0.43859E-01 ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 0.829315,  0.36934E-01 ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 0.808316,  0.35087E-01 ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 0.63903 , 0.21930E-01  ,7.485,  2.094,  0.55,	11.81,  63.54,  315,  12.00,  0.30000E+00},
								{ 3.26665 , 0.87700E-01  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.15450E+00},
								{ 2.00041 , 0.65760E-01  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.30000E+00},
								{ 1.70595 , 0.23920E-01  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.15430E+00},
								{ 1.41450 , 0.32880E-01  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.27300E+00},
								{ 1.29803 , 0.13850E-01  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.15430E+00},
								{ 1.15493 , 0.21925E-01  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.27270E+00},
								{ 1.08888 , 0.97400E-02  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.27270E+00},
								{ 0.86284 , 0.61200E-02  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.27270E+00},
								{ 0.79228 , 0.51588E-02  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.27270E+00},
								{ 0.79228 , 0.51600E-02  ,8.42 , 1.216,  0.18,  22.63,  72.6,  290,  9.0,  0.27270E+00},
								{ 3.13536 , 0.25970E-01  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.18080E+00},
								{ 1.92001 , 0.19480E-01  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.30000E+00},
								{ 1.63739 , 0.70800E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 1.35765 , 0.97400E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 1.24587 , 0.41000E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.18080E+00},
								{ 1.10852 , 0.64930E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 1.04512 , 0.28900E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 1.04512 , 0.28900E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 0.96000 , 0.48700E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 0.76044 , 0.15277E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 0.76044 , 0.15277E-02  ,2.18 , 0.0889,	0.0,  20.02,  28.09,  420,  6.36,  0.28000E+00},
								{ 1.97956 , 0.11361      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.30000E+00},
								{ 1.97956 , 0.11361      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.28000E+00},
								{ 1.32857 , 0.05117      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.28000E+00},
								{ 1.02290 , 0.091        ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.28000E+00},
								{ 1.14290 , 0.15147      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.28000E+00},
								{ 0.96363 , 0.10768      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.28000E+00},
								{ 0.98978 , 0.0284       ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.28000E+00},
								{ 1.79215 , 0.37245      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.30000E+00},
								{ 1.73285 , 0.26116      ,7.62579,  0.00422655,	0.002,  8.10926,  9.012,  1100,  7.62,  0.30000E+00},
								{ 3.35500 , 0.79500E+00  ,5.555,  0.0019,	0.0,  8.80,  12.01,  1050,  20.00,  0.30000E+00},
								{ 1.67750 , 0.18000E+00  ,5.555,  0.0019,	0.0,  8.80,  12.01,  1050,  20.00,  0.30000E+00},
								{ 1.11830 , 0.08833E+00  ,5.555,  0.0019,	0.0,  8.80,  12.01,  1050,  20.00,  0.30000E+00},
								{ 2.02660 , 0.34031E+00  ,11.43,  2.53,	0.4 , 11.75 , 55.85,  411,  10.67 , 0.30000E+00},
								{ 3.43500 , 0.11020E+00  ,1.79,  2.88,	0.55,  13.16,  48.0,  300,  12.00 , 0.30000E+00},
								{ 1.71750 , 0.13130E+00  ,1.79,  2.88,	0.55,  13.16,  48.0,  300,  12.00 , 0.30000E+00},
								{ 3.43500 , 0.55100E-01  ,1.79,  2.88,	0.55,  13.16,  48.0,  300,  12.00 , 0.30000E+00},
								{ 2.05929 , 0.36606      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 1.26105 , 0.27455      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 1.07543 , 0.09984      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.89170 , 0.13727      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.81828 , 0.0578       ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.72807 , 0.09152      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.68643 , 0.04067      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.68643 , 0.04067      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.63053 , 0.06864      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00},
								{ 0.63053 , 0.06864      ,5.55449  ,0.00194444,	0.0,  5.67213,  12.01,  1860,  3.00,  0.30000E+00}
	};

//////////////////////
// Testing function //
//////////////////////

	// Function to test what neutron status is throughout the simulation 
	void print_neutron_state(struct neutron_values* neutron){
		printf("Neutron state:\nki %g, %g, %g\nkf %g, %g, %g\nv %g, %g, %g\nr %g, %g, %g\nG %g, %g, %g\nki size %g, kf size %g, v size %g, wavelength %g, energy %g, Gr %g\n\n", 
			neutron->ki[0], neutron->ki[1], neutron->ki[2],
			neutron->kf[0], neutron->kf[1], neutron->kf[2],
			neutron->v[0], neutron->v[1], neutron->v[2],
			neutron->r[0], neutron->r[1], neutron->r[2],
			neutron->Gr[0], neutron->Gr[1], neutron->Gr[2],
			neutron->ki_size, neutron->kf_size, neutron->v_size, neutron->l, neutron->E, neutron->Gr_size
			);
	};

///////////////////////
// Neutron Functions //
///////////////////////
	// Set neutron values for calculation
	void set_neutron_values(
		struct neutron_values* neutron,
		double x, double y, double z,
		double vx, double vy, double vz){
			neutron->r[0] = x;
			neutron->r[1] = y;
			neutron->r[2] = z;
			neutron->v[0] = vx;
			neutron->v[1] = vy;
			neutron->v[2] = vz;
			neutron->v_size = 0;
			neutron->ki_size = 0;
			neutron->kf_size = 0;
			for (int i =0; i<3; i++){ 
				neutron->ki[i] = neutron->v[i]*V2K;
				neutron->ki_size += square(neutron->ki[i]);
				neutron->v_size += square(neutron->v[i]);
				neutron->kf_size += square(neutron->kf[i]);
			}        
			neutron->v_size = sqrt(neutron->v_size);
			neutron->ki_size = sqrt(neutron->ki_size);
			neutron->kf_size = sqrt(neutron->kf_size);
			neutron->l = 3956/neutron->v_size;// Wavelength in Angstrom
			neutron->E = square(neutron->v_size/437); // Energy in meV
	};

	// Find intersection points with monochromator.
	// It is an integer to check for edge effects.
	void calc_intersection(struct monochromator_values* monochromator, struct neutron_values* neutron){
		double translated_x;
		double inner_t0;
		double inner_t1;
		double outer_t0;
		double outer_t1;
		
		// Translate neutron to coordinates of center of radius of curvature
		if ((monochromator->type == 0 )||( monochromator->type == 2)){
			translated_x = neutron->r[0]-1000000000-monochromator->thickness/2;
		} else {
			translated_x = neutron->r[0]-monochromator->radius_horizontal-monochromator->thickness/2;
		}

		// Checks intersection for inner radius
		cylinder_intersect(&inner_t0,&inner_t1,
							translated_x,neutron->r[1],neutron->r[2],
							neutron->v[0],neutron->v[1],neutron->v[2],
							monochromator->radius_horizontal,
							monochromator->height);

		// Checks intersection for outer radius
		cylinder_intersect(&outer_t0,&outer_t1,
							translated_x,neutron->r[1],neutron->r[2],
							neutron->v[0],neutron->v[1],neutron->v[2],
							monochromator->radius_horizontal+monochromator->thickness,
							monochromator->height);
		
		// Defines which time to take. Will always scatter from the curvature furthest away,
		// so no "inverted monochromators". This also means the monochromators HAVE to be pointing towards
		// way the neutrons are coming. 
		if (neutron->direction == 1){
			if (neutron->scatter_count > 0){
				if (outer_t1 > 0){
					neutron->exit_time = outer_t0 > outer_t1 ? outer_t0 : outer_t1;
				} else {
					neutron->exit_time = 0;
				}
			} else {
				neutron->entry_time = inner_t0 > inner_t1 ? inner_t0 : inner_t1;
				neutron->exit_time = outer_t0 > outer_t1 ? outer_t0 : outer_t1;
			}
		}
		else {
			if ((inner_t0 > 0)){
				neutron->entry_time = inner_t0 < inner_t1 ? inner_t0 : inner_t1;
			} else {
				neutron->entry_time = 0;
			}
		}
	}

	// Calculate scattering vector G using Saroun eq. 3
	void calculate_G(struct monochromator_values* monochromator, struct neutron_values* neutron){
		neutron->Gr_size = 0;
		neutron->Grad_G_size = 0;
		double normalize = 0;
		// Curvature effect calculations
		for (int i=0; i<3; i++){
			neutron->Gr[i] = monochromator->G0[i];
			neutron->Gr_size += square(neutron->Gr[i]);
			for (int j=0; j<3; j++){
				neutron->Gr[i] += monochromator->lattice_spacing_gradient_field[i][j]*neutron->r[j];
			}
			neutron->Grad_G[i] = (neutron->Gr[i] - monochromator->G0[i])*neutron->direction;
			neutron->Grad_G_size += square(neutron->Grad_G[i]);
			neutron->Gr[i] *= neutron->direction;
		}
		// Set size for scattering vector
		neutron->Gr_size = sqrt(neutron->Gr_size);
		neutron->Grad_G_size = sqrt(neutron->Grad_G_size);


		// Define gamma vector, which has to be out of the ki-G plane
		vec_prod(neutron->perp_to_Gk[0],neutron->perp_to_Gk[1],neutron->perp_to_Gk[2],neutron->ki[0],neutron->ki[1],neutron->ki[2],neutron->Gr[0],neutron->Gr[1],neutron->Gr[2]);
		neutron->perp_size = sqrt(square(neutron->perp_to_Gk[0])+square(neutron->perp_to_Gk[1])+square(neutron->perp_to_Gk[2]));
		
		// Normal distribution of mosaicity
		if ((monochromator->type==2) || (monochromator->type==3)){
			neutron->gamma_0 = gauss_0(monochromator->mos);
		} else {
			neutron->gamma_0 = 0;
		}
		// Calculate Bragg angle from G and ki  
		neutron->Bragg_angle = PI - acos((neutron->ki[0]*monochromator->perp_G[0]+
									neutron->ki[1]*monochromator->perp_G[1]+
									neutron->ki[2]*monochromator->perp_G[2])/
									(neutron->ki_size));

		// Claculate Gamma vector
		for (int i=0;i<3;i++){
			neutron->perp_to_Gk[i] /= neutron->perp_size;
			neutron->Ggamma[i] = neutron->Gr_size*neutron->perp_to_Gk[i]*neutron->gamma_0;
			neutron->Gr[i] += neutron->Ggamma[i];	
			normalize += square(neutron->Gr[i]);					   
		}
		normalize = sqrt(normalize);
		neutron->Gr[0] *= neutron->Gr_size/normalize;
		neutron->Gr[1] *= neutron->Gr_size/normalize;
		neutron->Gr[2] *= neutron->Gr_size/normalize;
	}

	// Calculate epsilon zero and beta fron Saroun Eq. 4
	void calculate_epszero_and_beta(struct monochromator_values* monochromator, struct neutron_values* neutron){
		double G[3];
		double temp_eps0 = 0;
		double temp_beta = 0;
		G[0] = (neutron->ki[0]+monochromator->G0[0]*neutron->direction+neutron->Ggamma[0]);
		G[1] = (neutron->ki[1]+monochromator->G0[1]*neutron->direction+neutron->Ggamma[1]);
		G[2] = (neutron->ki[2]+monochromator->G0[2]*neutron->direction+neutron->Ggamma[2]);
		
		temp_eps0 = (square(G[1])+square(G[0])+square(G[2]) + square(neutron->Grad_G_size) - square(neutron->ki_size));
		double denominator = 2*neutron->Gr_size*neutron->ki_size*cos(neutron->Bragg_angle);
		neutron->eps_zero = -temp_eps0/denominator*neutron->direction;
		
		double Gk[3];
		for (int i=0;i<3;i++){
			Gk[i] = 0;
			for (int j=0;j<3;j++){
				Gk[i] += monochromator->lattice_spacing_gradient_field[i][j]*neutron->ki[j]*neutron->direction;
			}	
		}
		
		temp_beta = (G[0])*Gk[0]
					+ (G[1])*Gk[1]
					+ (G[2])*Gk[2];
		neutron->beta = -temp_beta/(denominator*neutron->ki_size/2);
	}

	// B0 and BT are values used for the Debye factor #TODO: Explain equation
	void calculate_B0_and_BT(struct monochromator_values *monochromator){
		double x;
		monochromator->B0 = 2872.556/monochromator->atomic_number
								/monochromator->debye_temperature;
		
		if (monochromator->temperature_mono>0.1) x = monochromator->debye_temperature/monochromator->temperature_mono;
		else x =monochromator->debye_temperature/0.1;
		double phi = calculate_phi_integral(x);

		monochromator->BT = 4 * monochromator->B0 * phi / square(x);
	}

	// #TODO: Find citation for this equation and check if it works
	// Calculate the kinematic reflectivity of material
	double calculate_kinematic_reflectivity(struct monochromator_values* monochromator, 
											struct neutron_values* neutron){
		double sine_of_bragg_angle = neutron->l/2/monochromator->lattice_spacing;
		double cosine_of_bragg_angle = sqrt(1-square(sine_of_bragg_angle));
		double extinction_length =  monochromator->lattice_spacing/neutron->l*sqrt(4*PI/monochromator->Maier_Leibnitz_reflectivity*100);
		
		// Kinenatic reflectivity = QML*DHKL*sin(theta_B)**2/PI/cos(theta_B)
		double kinematic_reflectivity = monochromator->Maier_Leibnitz_reflectivity;
		kinematic_reflectivity *= monochromator->lattice_spacing;
		kinematic_reflectivity *= square(sine_of_bragg_angle);
		kinematic_reflectivity *= 1/PI/cosine_of_bragg_angle;
		kinematic_reflectivity *= monochromator->Debye_Waller_factor;
		// Primary extinction factor, using the approximation in G.E Bacon and R.D. Lowde, Acta Cryst. (1948). 1, 303
		kinematic_reflectivity *= tanh(monochromator->domain_thickness/extinction_length)/monochromator->domain_thickness*extinction_length;
		return kinematic_reflectivity;
	}

	// Calculate time of flight and probability through crystal
	void calculate_tau_and_prob(struct monochromator_values* monochromator, struct neutron_values* neutron, int i){
		double Q = calculate_kinematic_reflectivity(monochromator, neutron);

		// First calculate probability of scattering
		if ((monochromator->type==2 )|| (monochromator->type==3)){
			neutron->P[i] = 1;
			double arg1 = (neutron->eps_zero+neutron->beta*neutron->max_tau*neutron->v_size)/monochromator->mos;
			neutron->P[i] = (1-exp(-(Q/(fabs(neutron->beta*neutron->ki_size)))*(normalCDF(arg1)-normalCDF(neutron->eps_zero/monochromator->mos))));

			// Second calculate time travelled through crystal
			neutron->tau[i] = normalCDF(neutron->eps_zero/monochromator->mos);
			neutron->tau[i]-= ((fabs(neutron->beta*neutron->ki_size))/Q)*log(1-rand01());
			neutron->tau[i] = inverse_of_normal_cdf(neutron->tau[i]);
			neutron->tau[i]*= monochromator->mos/(neutron->ki_size*neutron->beta);
			neutron->tau[i]-= neutron->eps_zero/(neutron->ki_size*neutron->beta);
			neutron->tau[i]*= neutron->ki_size/neutron->v_size;	
		} else {
			neutron->P[i] = 1- exp(-Q/(fabs(neutron->beta*neutron->ki_size)));
			neutron->tau[i] = -neutron->eps_zero/(neutron->ki_size*neutron->beta);
			neutron->tau[i] *= neutron->ki_size/neutron->v_size;
		}
	}

	// Update to final kf
	void calculate_kf(struct monochromator_values* monochromator, struct neutron_values* neutron, int j){

		// Calculate in G-ki plane mosaicity (Saroun eq. 4)
		double epsilon[3];
		double epsilon_size;
		double eps;
		eps = neutron->eps_zero+neutron->tau[j]*neutron->beta*neutron->v_size;
		// 	("%.13g, %.13g, %.13g, %.13g\n", neutron->eps_zero, neutron->beta, eps, neutron->P[j]);
		vec_prod(epsilon[0],epsilon[1],epsilon[2],neutron->perp_to_Gk[0],neutron->perp_to_Gk[1],neutron->perp_to_Gk[2],neutron->Gr[0],neutron->Gr[1],neutron->Gr[2]);
		epsilon_size = sqrt(square(epsilon[0])+square(epsilon[1])+square(epsilon[2]));
		neutron->kf_size = 0;

		// Calculate G
		double normalize = 0;
		for (int i=0;i<3;i++){
			neutron->Gr[i] = monochromator->G0[i];
			for (int k=0;k<3;k++){
				neutron->Gr[i] += monochromator->lattice_spacing_gradient_field[i][k]*neutron->r[k];
			}
			neutron->Gr[i] += neutron->direction*(neutron->Gr_size*eps*epsilon[i]/epsilon_size)+neutron->Ggamma[i];
			normalize += square(neutron->Gr[i]);
		}
		
		normalize = sqrt(normalize);
		for (int i=0;i<3;i++){
			neutron->Gr[i] *= neutron->Gr_size/normalize;
		}

		// Calculate final scattering vector (Saroun eq. 3)
		for (int i=0;i<3;i++){
			neutron->kf[i] = neutron->ki[i]+neutron->Gr[i];
			neutron->kf_size += square(neutron->kf[i]);
		}
		neutron->kf_size = sqrt(neutron->kf_size);

		// Noramlize for Bragg condition
		for (int i=0;i<3;i++){
			neutron->kf[i] *= neutron->ki_size/neutron->kf_size;
		}
	}

	// Reflect neutron
	void reflect_neutron(struct neutron_values* neutron, double* vx, double* vy, double* vz){
			*vx = (neutron->kf[0])*K2V;
			*vy = (neutron->kf[1])*K2V;
			*vz = (neutron->kf[2])*K2V;
	}

	// Attenuation coefficient calculation 
	// Cited from Andreas K. Freund 20. December 1982
	double calculate_attenuation_coefficient(struct monochromator_values* monochromator,
										struct neutron_values* neutron){
		double E = square(neutron->v_size)*VS2E; // Neutron energy in meV
		// Get factor for single phonon cross section
		
		double Bernoulli_sequence[31] = {1,-0.5,0.166667,0,-0.033333,0,0.0238095,0,-0.033333,
										0,0.0757576,0,-0.253114,0,1.16667,0,-7.09216,0,54.9712,
										0,-529.124,0,6192.12,0,-86580.3,0,1.42551717e6,0,-2.7298231e7,
										0,6.01580874e8};
		double x;
		if (monochromator->temperature_mono - 0.1 <= 0){
			x = monochromator->debye_temperature/0.1;
		}
		else{
			x = monochromator->debye_temperature/monochromator->temperature_mono;
		} 
		double R, Ifact, Xn;
		if (x<6){
			R = 0;
			Ifact = 1;
			Xn = 1/x;
			for (int i=0; i<30; i++){
				R += Bernoulli_sequence[i]*Xn/Ifact/(i + 2.5);
				Xn *= x;
				Ifact *= i + 1;
			}
		}
		else R = 3.3/sqrt(x*x*x*x*x*x*x);

		// Define boltzmann_constant in units of (meV/K)
		double boltzmann_constant = 0.08617333262;
		double DWMF =  1-exp(-(monochromator->B0+monochromator->BT)*monochromator->Constant_from_Freund_paper*E/1000);
		// Set the cross sections, as written in freunds paper
		monochromator->nuclear_capture_absorption = monochromator->incoherent_scattering_cross_section
													+monochromator->absorption_for_1AA_Neutrons*neutron->l;

		monochromator->multiple_phonon_absorption = monochromator->bound_atom_scattering_cross_section
											*square(monochromator->atomic_number/(monochromator->atomic_number + 1))
											*DWMF;

		monochromator->single_phonon_absorption = 3*monochromator->bound_atom_scattering_cross_section/monochromator->atomic_number
											* sqrt(boltzmann_constant * monochromator->debye_temperature/E) * R;

		double attenuation_coefficient =  (monochromator->nuclear_capture_absorption
										+ monochromator->single_phonon_absorption
										+ monochromator->multiple_phonon_absorption)
										/monochromator->volume * 100; // *100 to change to per mm?
		return attenuation_coefficient;
	}

	void prop_half(struct monochromator_values* monochromator, struct neutron_values* neutron){
		// Calculate intersection with monochromator
		neutron->r[0] += neutron->v[0]*(neutron->exit_time-neutron->entry_time)/2;
		neutron->r[1] += neutron->v[1]*(neutron->exit_time-neutron->entry_time)/2;
		neutron->r[2] += neutron->v[2]*(neutron->exit_time-neutron->entry_time)/2;
	}


/* Shared user declarations for all components types 'Monitor_nD'. */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Modified by: TW, Nov 2020: introduced user doubles
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H

#define MONITOR_ND_LIB_H "$Revision$"
#define MONnD_COORD_NMAX  30  /* max number of variables to record */

  typedef struct MonitornD_Defines
  {
    int COORD_NONE  ;
    int COORD_X     ;
    int COORD_Y     ;
    int COORD_Z     ;
    int COORD_RADIUS; 
    int COORD_VX    ;
    int COORD_VY    ;
    int COORD_VZ    ;
    int COORD_V     ;
    int COORD_T     ;
    int COORD_P     ;
    int COORD_SX    ;
    int COORD_SY    ;
    int COORD_SZ    ;
    int COORD_KX    ;
    int COORD_KY    ;
    int COORD_KZ    ;
    int COORD_K     ;
    int COORD_ENERGY;
    int COORD_LAMBDA;
    int COORD_KXY   ;
    int COORD_KYZ   ;
    int COORD_KXZ   ;
    int COORD_VXY   ;
    int COORD_VYZ   ;
    int COORD_VXZ   ;
    int COORD_HDIV  ;
    int COORD_VDIV  ;
    int COORD_ANGLE ;
    int COORD_NCOUNT;
    int COORD_THETA ;
    int COORD_PHI   ;
    int COORD_USER1 ;
    int COORD_USER2 ;
    int COORD_USER3 ;
    int COORD_USERDOUBLE0 ;
    int COORD_USERDOUBLE1 ;
    int COORD_USERDOUBLE2 ;
    int COORD_USERDOUBLE3 ;
    int COORD_USERDOUBLE4 ;
    int COORD_USERDOUBLE5 ;
    int COORD_USERDOUBLE6 ;
    int COORD_USERDOUBLE7 ;
    int COORD_USERDOUBLE8 ;
    int COORD_USERDOUBLE9 ;
    int COORD_USERDOUBLE10 ;
    int COORD_USERDOUBLE11 ;
    int COORD_USERDOUBLE12 ;
    int COORD_USERDOUBLE13 ;
    int COORD_USERDOUBLE14 ;
    int COORD_USERDOUBLE15 ;
    int COORD_XY    ;
    int COORD_XZ    ;
    int COORD_YZ    ;
    int COORD_PIXELID;

    /* token modifiers */
    int COORD_VAR   ; /* next token should be a variable or normal option */
    int COORD_MIN   ; /* next token is a min value */
    int COORD_MAX   ; /* next token is a max value */
    int COORD_DIM   ; /* next token is a bin value */
    int COORD_FIL   ; /* next token is a filename */
    int COORD_EVNT  ; /* next token is a buffer size value */
    int COORD_3HE   ; /* next token is a 3He pressure value */
    int COORD_LOG   ; /* next variable will be in log scale */
    int COORD_ABS   ; /* next variable will be in abs scale */
    int COORD_SIGNAL; /* next variable will be the signal var */
    int COORD_AUTO  ; /* set auto limits */

    char TOKEN_DEL[32]; /* token separators */

    char SHAPE_SQUARE; /* shape of the monitor */
    char SHAPE_DISK  ;
    char SHAPE_SPHERE;
    char SHAPE_CYLIND;
    char SHAPE_BANANA; /* cylinder without top/bottom, on restricted angular area */
    char SHAPE_BOX   ;
    char SHAPE_PREVIOUS;
    char SHAPE_OFF;

  } MonitornD_Defines_type;

  typedef struct MonitornD_Variables
  {
    double area;
    double Sphere_Radius     ;
    double Cylinder_Height   ;
    char   Flag_With_Borders ;   /* 2 means xy borders too */
    char   Flag_List         ;   /* 1 store 1 buffer, 2 is list all, 3 list all+append */
    char   Flag_Multiple     ;   /* 1 when n1D, 0 for 2D */
    char   Flag_Verbose      ;
    int    Flag_Shape        ;
    char   Flag_Auto_Limits  ;   /* get limits from first Buffer */
    char   Flag_Absorb       ;   /* monitor is also a slit */
    char   Flag_per_cm2      ;   /* flux is per cm2 */
    char   Flag_log          ;   /* log10 of the flux */
    char   Flag_parallel     ;   /* set neutron state back after detection (parallel components) */
    char   Flag_Binary_List  ;
    char   Flag_capture      ;   /* lambda monitor with lambda/lambda(2200m/s = 1.7985 Angs) weightening */
    int    Flag_signal       ;   /* 0:monitor p, else monitor a mean value */
    int    Flag_mantid       ;   /* 0:normal monitor, else do mantid-event specifics */
    int    Flag_OFF          ;   /* Flag to indicate external geometry from OFF file */
    long long OFF_polyidx;   /* When intersection is done externally by off_intersect, this gives the 
				    polygon number, i.e. pixel index */

    unsigned long Coord_Number      ;   /* total number of variables to monitor, plus intensity (0) */
    unsigned long Coord_NumberNoPixel;  /* same but without counting PixelID */
    unsigned long Buffer_Block      ;   /* Buffer size for list or auto limits */
    long long Neutron_Counter   ;   /* event counter, simulation total counts is mcget_ncount() */
    unsigned long Buffer_Counter    ;   /* index in Buffer size (for realloc) */
    unsigned long Buffer_Size       ;
    int    Coord_Type[MONnD_COORD_NMAX];      /* type of variable */
    char   Coord_Label[MONnD_COORD_NMAX][30]; /* label of variable */
    char   Coord_Var[MONnD_COORD_NMAX][30];   /* short id of variable */
    long   Coord_Bin[MONnD_COORD_NMAX];       /* bins of variable array */
    long   Coord_BinProd[MONnD_COORD_NMAX];   /* product of bins of variable array */
    double Coord_Min[MONnD_COORD_NMAX];
    double Coord_Max[MONnD_COORD_NMAX];
    char   Monitor_Label[MONnD_COORD_NMAX*30];/* Label for monitor */
    char   Mon_File[128];                     /* output file name */

    /* these don't seem to be used anymore as they are superseded by _particle
    double cx, cy, cz;
    double cvx, cvy, cvz;
    double ckx, cky, ckz;
    double csx, csy, csz;
    double cEx, cEy, cEz;
    double cs1, cs2, ct, cphi, cp; */

    double He3_pressure;
    char   Flag_UsePreMonitor    ;   /* use a previously stored neutron parameter set */
    char   UserName1[128];
    char   UserName2[128];
    char   UserName3[128];
    char   UserVariable1[128];
    char   UserVariable2[128];
    char   UserVariable3[128];
    double UserDoubles[16];
    char   option[CHAR_BUF_LENGTH];

    long long int Nsum;
    double psum, p2sum;
    double **Mon2D_N;
    double **Mon2D_p;
    double **Mon2D_p2;
    double *Mon2D_Buffer;
    unsigned long PixelID;

    double mxmin,mxmax,mymin,mymax,mzmin,mzmax;
    double mean_dx, mean_dy, min_x, min_y, max_x, max_y, mean_p;

    char   compcurname[128];
    Coords compcurpos;

  } MonitornD_Variables_type;

/* monitor_nd-lib function prototypes */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *, MonitornD_Variables_type *, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, MCNUM, int);
#pragma acc routine
int Monitor_nD_Trace(MonitornD_Defines_type *, MonitornD_Variables_type *, _class_particle* _particle);
MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_Finally(MonitornD_Defines_type *, MonitornD_Variables_type *);
void Monitor_nD_McDisplay(MonitornD_Defines_type *, MonitornD_Variables_type *);

#endif

/* end of monitor_nd-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/monitor_nd-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Modified by: TW, Nov 2020: introduced user doubles
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by the monitor_nd related components
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "monitor_nd-lib"
*
*******************************************************************************/

#ifndef MONITOR_ND_LIB_H
#error McStas : please import this library with %include "monitor_nd-lib"
#endif

/* ========================================================================= */
/* Monitor_nD_Init: this routine is used to parse options                    */
/* ========================================================================= */

void Monitor_nD_Init(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars,
  MCNUM xwidth,
  MCNUM yheight,
  MCNUM zdepth,
  MCNUM xmin,
  MCNUM xmax,
  MCNUM ymin,
  MCNUM ymax,
  MCNUM zmin,
  MCNUM zmax,
  int offflag)
  {
    long carg = 1;
    char *option_copy, *token;
    char Flag_New_token = 1;
    char Flag_End       = 1;
    char Flag_All       = 0;
    char Flag_No        = 0;
    char Flag_abs       = 0;
    int  Flag_auto      = 0;  /* -1: all, 1: the current variable */
    int  Set_Vars_Coord_Type;
    char Set_Vars_Coord_Label[64];
    char Set_Vars_Coord_Var[64];
    char Short_Label[MONnD_COORD_NMAX][64];
    int  Set_Coord_Mode;
    long i=0, j=0;
    double lmin, lmax, XY=0;
    long t;


    t = (long)time(NULL);

/* initialize DEFS */
/* Variables to monitor */
    DEFS->COORD_NONE   =0;
    DEFS->COORD_X      =1;
    DEFS->COORD_Y      =2;
    DEFS->COORD_Z      =3;
    DEFS->COORD_RADIUS =19;
    DEFS->COORD_VX     =4;
    DEFS->COORD_VY     =5;
    DEFS->COORD_VZ     =6;
    DEFS->COORD_V      =16;
    DEFS->COORD_T      =7;
    DEFS->COORD_P      =8;
    DEFS->COORD_SX     =9;
    DEFS->COORD_SY     =10;
    DEFS->COORD_SZ     =11;
    DEFS->COORD_KX     =12;
    DEFS->COORD_KY     =13;
    DEFS->COORD_KZ     =14;
    DEFS->COORD_K      =15;
    DEFS->COORD_ENERGY =17;
    DEFS->COORD_LAMBDA =18;
    DEFS->COORD_HDIV   =20;
    DEFS->COORD_VDIV   =21;
    DEFS->COORD_ANGLE  =22;
    DEFS->COORD_NCOUNT =23;
    DEFS->COORD_THETA  =24;
    DEFS->COORD_PHI    =25;
    DEFS->COORD_USER1  =26;
    DEFS->COORD_USER2  =27;
    DEFS->COORD_USER3  =28;
    DEFS->COORD_USERDOUBLE0=39;
    DEFS->COORD_USERDOUBLE1=40;
    DEFS->COORD_USERDOUBLE2=41;
    DEFS->COORD_USERDOUBLE3=42;
    DEFS->COORD_USERDOUBLE4=43;
    DEFS->COORD_USERDOUBLE5=44;
    DEFS->COORD_USERDOUBLE6=45;
    DEFS->COORD_USERDOUBLE7=46;
    DEFS->COORD_USERDOUBLE8=47;
    DEFS->COORD_USERDOUBLE9=48;
    DEFS->COORD_USERDOUBLE10=49;
    DEFS->COORD_USERDOUBLE11=50;
    DEFS->COORD_USERDOUBLE12=51;
    DEFS->COORD_USERDOUBLE13=52;
    DEFS->COORD_USERDOUBLE14=53;
    DEFS->COORD_USERDOUBLE15=54;
    DEFS->COORD_XY     =37;
    DEFS->COORD_YZ     =31;
    DEFS->COORD_XZ     =32;
    DEFS->COORD_VXY    =30;
    DEFS->COORD_VYZ    =34;
    DEFS->COORD_VXZ    =36;
    DEFS->COORD_KXY    =29;
    DEFS->COORD_KYZ    =33;
    DEFS->COORD_KXZ    =35;
    DEFS->COORD_PIXELID=38;

/* token modifiers */
    DEFS->COORD_VAR    =0;    /* next token should be a variable or normal option */
    DEFS->COORD_MIN    =1;    /* next token is a min value */
    DEFS->COORD_MAX    =2;    /* next token is a max value */
    DEFS->COORD_DIM    =3;    /* next token is a bin value */
    DEFS->COORD_FIL    =4;    /* next token is a filename */
    DEFS->COORD_EVNT   =5;    /* next token is a buffer size value */
    DEFS->COORD_3HE    =6;    /* next token is a 3He pressure value */
    DEFS->COORD_LOG    =64;   /* next variable will be in log scale */
    DEFS->COORD_ABS    =128;  /* next variable will be in abs scale */
    DEFS->COORD_SIGNAL =256;  /* next variable will be the signal var */
    DEFS->COORD_AUTO   =512;  /* set auto limits */

    strcpy(DEFS->TOKEN_DEL, " =,;[](){}:");  /* token separators */

    DEFS->SHAPE_SQUARE =0;    /* shape of the monitor */
    DEFS->SHAPE_DISK   =1;
    DEFS->SHAPE_SPHERE =2;
    DEFS->SHAPE_CYLIND =3;
    DEFS->SHAPE_BANANA =4;
    DEFS->SHAPE_BOX    =5;
    DEFS->SHAPE_PREVIOUS=6;
    DEFS->SHAPE_OFF=7;

    Vars->Sphere_Radius     = 0;
    Vars->Cylinder_Height   = 0;
    Vars->Flag_With_Borders = 0;   /* 2 means xy borders too */
    Vars->Flag_List         = 0;   /* 1=store 1 buffer, 2=list all, 3=re-use buffer */
    Vars->Flag_Multiple     = 0;   /* 1 when n1D, 0 for 2D */
    Vars->Flag_Verbose      = 0;
    Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    Vars->Flag_Auto_Limits  = 0;   /* get limits from first Buffer */
    Vars->Flag_Absorb       = 0;   /* monitor is also a slit */
    Vars->Flag_per_cm2      = 0;   /* flux is per cm2 */
    Vars->Flag_log          = 0;   /* log10 of the flux */
    Vars->Flag_parallel     = 0;   /* set neutron state back after detection (parallel components) */
    Vars->Flag_Binary_List  = 0;   /* save list as a binary file (smaller) */
    Vars->Coord_Number      = 0;   /* total number of variables to monitor, plus intensity (0) */
    Vars->Coord_NumberNoPixel=0;   /* same but without counting PixelID */

    Vars->Buffer_Block      = MONND_BUFSIZ;     /* Buffer size for list or auto limits */	
    Vars->Neutron_Counter   = -1;   /* event counter, simulation total counts is mcget_ncount() */
    Vars->Buffer_Counter    = 0;   /* index in Buffer size (for realloc) */
    Vars->Buffer_Size       = 0;
    Vars->He3_pressure      = 0;
    Vars->Flag_capture      = 0;
    Vars->Flag_signal       = DEFS->COORD_P;
    Vars->Flag_mantid       = 0;
    Vars->Flag_OFF          = offflag;
    Vars->OFF_polyidx       = -1;
    Vars->mean_dx=Vars->mean_dy=0;
    Vars->min_x = Vars->max_x  =0;
    Vars->min_y = Vars->max_y  =0;

    Set_Vars_Coord_Type = DEFS->COORD_NONE;
    Set_Coord_Mode = DEFS->COORD_VAR;

    /* handle size parameters */
    /* normal use is with xwidth, yheight, zdepth */
    /* if xmin,xmax,ymin,ymax,zmin,zmax are non 0, use them */
    if (fabs(xmin-xmax) == 0)
      { Vars->mxmin = -fabs(xwidth)/2; Vars->mxmax = fabs(xwidth)/2; }
    else
      { if (xmin < xmax) {Vars->mxmin = xmin; Vars->mxmax = xmax;}
        else {Vars->mxmin = xmax; Vars->mxmax = xmin;}
      }
    if (fabs(ymin-ymax) == 0)
      { Vars->mymin = -fabs(yheight)/2; Vars->mymax = fabs(yheight)/2; }
    else
      { if (ymin < ymax) {Vars->mymin = ymin; Vars->mymax = ymax;}
        else {Vars->mymin = ymax; Vars->mymax = ymin;}
      }
    if (fabs(zmin-zmax) == 0)
      { Vars->mzmin = -fabs(zdepth)/2; Vars->mzmax = fabs(zdepth)/2; }
    else
      { if (zmin < zmax) {Vars->mzmin = zmin; Vars->mzmax = zmax; }
        else {Vars->mzmin = zmax; Vars->mzmax = zmin; }
      }

    if (fabs(Vars->mzmax-Vars->mzmin) == 0)
      Vars->Flag_Shape        = DEFS->SHAPE_SQUARE;
    else
      Vars->Flag_Shape        = DEFS->SHAPE_BOX;

    if (Vars->Flag_OFF) {
      Vars->Flag_Shape        = DEFS->SHAPE_OFF;
    }
    
    /* parse option string */

    option_copy = (char*)malloc(strlen(Vars->option)+1);
    if (option_copy == NULL)
    {
      fprintf(stderr,"Monitor_nD: %s cannot allocate 'options' copy (%li). Fatal.\n", Vars->compcurname, (long)strlen(Vars->option));
      exit(-1);
    }

    if (strlen(Vars->option))
    {
      Flag_End = 0;
      strcpy(option_copy, Vars->option);
    }

    if (strstr(Vars->option, "cm2") || strstr(Vars->option, "cm^2")) Vars->Flag_per_cm2 = 1;

    if (strstr(Vars->option, "binary") || strstr(Vars->option, "float"))
      Vars->Flag_Binary_List  = 1;
    if (strstr(Vars->option, "double"))
      Vars->Flag_Binary_List  = 2;

    strcpy(Vars->Coord_Label[0],"Intensity");
    strncpy(Vars->Coord_Var[0],"p",30);
    Vars->Coord_Type[0] = DEFS->COORD_P;
    Vars->Coord_Bin[0] = 1;
    Vars->Coord_Min[0] = 0;
    Vars->Coord_Max[0] = FLT_MAX;

    /* default file name is comp_name+dateID */
    sprintf(Vars->Mon_File, "%s_%li", Vars->compcurname, t);

    carg = 1;
    while((Flag_End == 0) && (carg < 128))
    {

      if (Flag_New_token) /* retain previous token or get a new one */
      {
        if (carg == 1) token=(char *)strtok(option_copy,DEFS->TOKEN_DEL);
        else token=(char *)strtok(NULL,DEFS->TOKEN_DEL);
        if (token == NULL) Flag_End=1;
      }
      Flag_New_token = 1;
      if ((token != NULL) && (strlen(token) != 0))
      {
        char iskeyword=0; /* left at 0 when variables are processed, 1 for modifiers */
        int  old_Mode;
        /* change token to lower case */
        for (i=0; i<strlen(token); i++) token[i]=tolower(token[i]);
        /* first handle option values from preceeding keyword token detected */
        old_Mode = Set_Coord_Mode;
        if (Set_Coord_Mode == DEFS->COORD_MAX)  /* max=%i */
        {
          if (!Flag_All)
            Vars->Coord_Max[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Max[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_MIN)  /* min=%i */
        {
          if (!Flag_All)
            Vars->Coord_Min[Vars->Coord_Number] = atof(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Min[i++] = atof(token));
          Set_Coord_Mode = DEFS->COORD_MAX;
        }
        if (Set_Coord_Mode == DEFS->COORD_DIM)  /* bins=%i */
        {
          if (!Flag_All)
            Vars->Coord_Bin[Vars->Coord_Number] = atoi(token);
          else
            for (i = 0; i <= Vars->Coord_Number; Vars->Coord_Bin[i++] = atoi(token));
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_FIL)  /* file=%s */
        {
          if (!Flag_No) strncpy(Vars->Mon_File,token,128);
          else { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1;}
          Set_Coord_Mode = DEFS->COORD_VAR;
        }
        if (Set_Coord_Mode == DEFS->COORD_EVNT) /* list=%i */
        {
          if (!strcmp(token, "all") || Flag_All) Vars->Flag_List = 2;
          else { i = (long)ceil(atof(token)); if (i) Vars->Buffer_Block = i;
            Vars->Flag_List = 1; }
          Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }
        if (Set_Coord_Mode == DEFS->COORD_3HE)  /* pressure=%g */
        {
            Vars->He3_pressure = atof(token);
            Set_Coord_Mode = DEFS->COORD_VAR; Flag_All = 0;
        }

        /* now look for general option keywords */
        if (!strcmp(token, "borders"))  {Vars->Flag_With_Borders = 1; iskeyword=1; }
        if (!strcmp(token, "verbose"))  {Vars->Flag_Verbose      = 1; iskeyword=1; }
        if (!strcmp(token, "log"))      {Vars->Flag_log          = 1; iskeyword=1; }
        if (!strcmp(token, "abs"))      {Flag_abs                = 1; iskeyword=1; }
        if (!strcmp(token, "multiple")) {Vars->Flag_Multiple     = 1; iskeyword=1; }
        if (!strcmp(token, "list") || !strcmp(token, "events")) {
          Vars->Flag_List = 1; Set_Coord_Mode = DEFS->COORD_EVNT;  }
        if (!strcmp(token, "limits") || !strcmp(token, "min"))
          Set_Coord_Mode = DEFS->COORD_MIN;
        if (!strcmp(token, "slit") || !strcmp(token, "absorb")) {
          Vars->Flag_Absorb = 1;  iskeyword=1; }
        if (!strcmp(token, "max"))  Set_Coord_Mode = DEFS->COORD_MAX;
        if (!strcmp(token, "bins") || !strcmp(token, "dim")) Set_Coord_Mode = DEFS->COORD_DIM;
        if (!strcmp(token, "file") || !strcmp(token, "filename")) {
          Set_Coord_Mode = DEFS->COORD_FIL;
          if (Flag_No) { strcpy(Vars->Mon_File,""); Vars->Coord_Number = 0; Flag_End = 1; }
        }
        if (!strcmp(token, "unactivate")) {
          Flag_End = 1; Vars->Coord_Number = 0; iskeyword=1; }
        if (!strcmp(token, "all"))    { Flag_All = 1;  iskeyword=1; }
        if (!strcmp(token, "sphere")) { Vars->Flag_Shape = DEFS->SHAPE_SPHERE; iskeyword=1; }
        if (!strcmp(token, "cylinder")) { Vars->Flag_Shape = DEFS->SHAPE_CYLIND; iskeyword=1; }
        if (!strcmp(token, "banana")) { Vars->Flag_Shape = DEFS->SHAPE_BANANA; iskeyword=1; }
        if (!strcmp(token, "square")) { Vars->Flag_Shape = DEFS->SHAPE_SQUARE; iskeyword=1; }
        if (!strcmp(token, "disk"))   { Vars->Flag_Shape = DEFS->SHAPE_DISK; iskeyword=1; }
        if (!strcmp(token, "box"))     { Vars->Flag_Shape = DEFS->SHAPE_BOX; iskeyword=1; }
        if (!strcmp(token, "previous")) { Vars->Flag_Shape = DEFS->SHAPE_PREVIOUS; iskeyword=1; }
        if (!strcmp(token, "parallel")){ Vars->Flag_parallel = 1; iskeyword=1; }
        if (!strcmp(token, "capture")) { Vars->Flag_capture = 1; iskeyword=1; }
        if (!strcmp(token, "auto")) { 
        #ifndef OPENACC
        if (Flag_auto != -1) {
	    Vars->Flag_Auto_Limits = 1;
	    if (Flag_All) Flag_auto = -1;
	    else          Flag_auto = 1;
	    iskeyword=1; Flag_All=0; 
	  }
        #endif
	}
        if (!strcmp(token, "premonitor")) {
          Vars->Flag_UsePreMonitor = 1; iskeyword=1; }
        if (!strcmp(token, "3He_pressure") || !strcmp(token, "pressure")) {
          Vars->He3_pressure = 3; iskeyword=1; }
        if (!strcmp(token, "no") || !strcmp(token, "not")) { Flag_No = 1;  iskeyword=1; }
        if (!strcmp(token, "signal")) Set_Coord_Mode = DEFS->COORD_SIGNAL;
        if (!strcmp(token, "mantid")) { Vars->Flag_mantid = 1; iskeyword=1; }

        /* Mode has changed: this was a keyword or value  ? */
        if (Set_Coord_Mode != old_Mode) iskeyword=1;

        /* now look for variable names to monitor */
        Set_Vars_Coord_Type = DEFS->COORD_NONE; lmin = 0; lmax = 0;

        if (!strcmp(token, "x"))
          { Set_Vars_Coord_Type = DEFS->COORD_X; strcpy(Set_Vars_Coord_Label,"x [m]"); strcpy(Set_Vars_Coord_Var,"x");
          lmin = Vars->mxmin; lmax = Vars->mxmax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mxmin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mxmax;}
        if (!strcmp(token, "y"))
          { Set_Vars_Coord_Type = DEFS->COORD_Y; strcpy(Set_Vars_Coord_Label,"y [m]"); strcpy(Set_Vars_Coord_Var,"y");
          lmin = Vars->mymin; lmax = Vars->mymax;
          Vars->Coord_Min[Vars->Coord_Number+1] = Vars->mymin;
          Vars->Coord_Max[Vars->Coord_Number+1] = Vars->mymax;}
        if (!strcmp(token, "z"))
          { Set_Vars_Coord_Type = DEFS->COORD_Z; strcpy(Set_Vars_Coord_Label,"z [m]"); strcpy(Set_Vars_Coord_Var,"z"); lmin = Vars->mzmin; lmax = Vars->mzmax; }
        if (!strcmp(token, "k") || !strcmp(token, "wavevector"))
          { Set_Vars_Coord_Type = DEFS->COORD_K; strcpy(Set_Vars_Coord_Label,"|k| [Angs-1]"); strcpy(Set_Vars_Coord_Var,"k"); lmin = 0; lmax = 10; }
        if (!strcmp(token, "v"))
          { Set_Vars_Coord_Type = DEFS->COORD_V; strcpy(Set_Vars_Coord_Label,"Velocity [m/s]"); strcpy(Set_Vars_Coord_Var,"v"); lmin = 0; lmax = 10000; }
        if (!strcmp(token, "t") || !strcmp(token, "time") || !strcmp(token, "tof"))
          { Set_Vars_Coord_Type = DEFS->COORD_T; strcpy(Set_Vars_Coord_Label,"TOF [s]"); strcpy(Set_Vars_Coord_Var,"t"); lmin = 0; lmax = 1.0; }
        if ((!strcmp(token, "p") || !strcmp(token, "i") || !strcmp(token, "intensity") || !strcmp(token, "flux")))
          { Set_Vars_Coord_Type = DEFS->COORD_P;
            strcpy(Set_Vars_Coord_Label,"Intensity");
            strncat(Set_Vars_Coord_Label, " [n/s", 30);
            if (Vars->Flag_per_cm2) strncat(Set_Vars_Coord_Label, "/cm2", 30);
            if (XY > 1 && Vars->Coord_Number)
              strncat(Set_Vars_Coord_Label, "/bin", 30);
            strncat(Set_Vars_Coord_Label, "]", 30);
            strcpy(Set_Vars_Coord_Var,"I");
            lmin = 0; lmax = FLT_MAX;
            if (Flag_auto>0) Flag_auto=0;
          }

        if (!strcmp(token, "vx"))
          { Set_Vars_Coord_Type = DEFS->COORD_VX; strcpy(Set_Vars_Coord_Label,"vx [m/s]"); strcpy(Set_Vars_Coord_Var,"vx"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VY; strcpy(Set_Vars_Coord_Label,"vy [m/s]"); strcpy(Set_Vars_Coord_Var,"vy"); lmin = -1000; lmax = 1000; }
        if (!strcmp(token, "vz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VZ; strcpy(Set_Vars_Coord_Label,"vz [m/s]"); strcpy(Set_Vars_Coord_Var,"vz"); lmin = -10000; lmax = 10000; }
        if (!strcmp(token, "kx"))
          { Set_Vars_Coord_Type = DEFS->COORD_KX; strcpy(Set_Vars_Coord_Label,"kx [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "ky"))
          { Set_Vars_Coord_Type = DEFS->COORD_KY; strcpy(Set_Vars_Coord_Label,"ky [Angs-1]"); strcpy(Set_Vars_Coord_Var,"ky"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "kz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KZ; strcpy(Set_Vars_Coord_Label,"kz [Angs-1]"); strcpy(Set_Vars_Coord_Var,"kz"); lmin = -10; lmax = 10; }
        if (!strcmp(token, "sx"))
          { Set_Vars_Coord_Type = DEFS->COORD_SX; strcpy(Set_Vars_Coord_Label,"sx [1]"); strcpy(Set_Vars_Coord_Var,"sx"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sy"))
          { Set_Vars_Coord_Type = DEFS->COORD_SY; strcpy(Set_Vars_Coord_Label,"sy [1]"); strcpy(Set_Vars_Coord_Var,"sy"); lmin = -1; lmax = 1; }
        if (!strcmp(token, "sz"))
          { Set_Vars_Coord_Type = DEFS->COORD_SZ; strcpy(Set_Vars_Coord_Label,"sz [1]"); strcpy(Set_Vars_Coord_Var,"sz"); lmin = -1; lmax = 1; }

        if (!strcmp(token, "energy") || !strcmp(token, "omega") || !strcmp(token, "e"))
          { Set_Vars_Coord_Type = DEFS->COORD_ENERGY; strcpy(Set_Vars_Coord_Label,"Energy [meV]"); strcpy(Set_Vars_Coord_Var,"E"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "lambda") || !strcmp(token, "wavelength") || !strcmp(token, "l"))
          { Set_Vars_Coord_Type = DEFS->COORD_LAMBDA; strcpy(Set_Vars_Coord_Label,"Wavelength [Angs]"); strcpy(Set_Vars_Coord_Var,"L"); lmin = 0; lmax = 100; }
        if (!strcmp(token, "radius") || !strcmp(token, "r"))
          { Set_Vars_Coord_Type = DEFS->COORD_RADIUS; strcpy(Set_Vars_Coord_Label,"Radius [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xy"))
          { Set_Vars_Coord_Type = DEFS->COORD_XY; strcpy(Set_Vars_Coord_Label,"Radius (xy) [m]"); strcpy(Set_Vars_Coord_Var,"xy"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "yz"))
          { Set_Vars_Coord_Type = DEFS->COORD_YZ; strcpy(Set_Vars_Coord_Label,"Radius (yz) [m]"); strcpy(Set_Vars_Coord_Var,"yz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "xz"))
          { Set_Vars_Coord_Type = DEFS->COORD_XZ; strcpy(Set_Vars_Coord_Label,"Radius (xz) [m]"); strcpy(Set_Vars_Coord_Var,"xz"); lmin = 0; lmax = xmax; }
        if (!strcmp(token, "vxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXY; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xy) [m]"); strcpy(Set_Vars_Coord_Var,"Vxy"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxy"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXY; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xy) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxy"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VYZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (yz) [m]"); strcpy(Set_Vars_Coord_Var,"Vyz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kyz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KYZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (yz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kyz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "vxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_VXZ; strcpy(Set_Vars_Coord_Label,"Radial Velocity (xz) [m]"); strcpy(Set_Vars_Coord_Var,"Vxz"); lmin = 0; lmax = 2000; }
        if (!strcmp(token, "kxz"))
          { Set_Vars_Coord_Type = DEFS->COORD_KXZ; strcpy(Set_Vars_Coord_Label,"Radial Wavevector (xz) [Angs-1]"); strcpy(Set_Vars_Coord_Var,"Kxz"); lmin = 0; lmax = 2; }
        if (!strcmp(token, "angle") || !strcmp(token, "a"))
          { Set_Vars_Coord_Type = DEFS->COORD_ANGLE; strcpy(Set_Vars_Coord_Label,"Angle [deg]"); strcpy(Set_Vars_Coord_Var,"A"); lmin = -50; lmax = 50; }
        if (!strcmp(token, "hdiv")|| !strcmp(token, "divergence") || !strcmp(token, "xdiv") || !strcmp(token, "hd") || !strcmp(token, "dx"))
          { Set_Vars_Coord_Type = DEFS->COORD_HDIV; strcpy(Set_Vars_Coord_Label,"Hor. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"hd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "vdiv") || !strcmp(token, "ydiv") || !strcmp(token, "vd") || !strcmp(token, "dy"))
          { Set_Vars_Coord_Type = DEFS->COORD_VDIV; strcpy(Set_Vars_Coord_Label,"Vert. Divergence [deg]"); strcpy(Set_Vars_Coord_Var,"vd"); lmin = -5; lmax = 5; }
        if (!strcmp(token, "theta") || !strcmp(token, "longitude") || !strcmp(token, "th"))
          { Set_Vars_Coord_Type = DEFS->COORD_THETA; strcpy(Set_Vars_Coord_Label,"Longitude [deg]"); strcpy(Set_Vars_Coord_Var,"th"); lmin = -180; lmax = 180; }
        if (!strcmp(token, "phi") || !strcmp(token, "latitude") || !strcmp(token, "ph"))
          { Set_Vars_Coord_Type = DEFS->COORD_PHI; strcpy(Set_Vars_Coord_Label,"Latitude [deg]"); strcpy(Set_Vars_Coord_Var,"ph"); lmin = -90; lmax = 90; }
        if (!strcmp(token, "ncounts") || !strcmp(token, "n") || !strcmp(token, "neutron"))
          { Set_Vars_Coord_Type = DEFS->COORD_NCOUNT; strcpy(Set_Vars_Coord_Label,"Neutron ID [1]"); strcpy(Set_Vars_Coord_Var,"n"); lmin = 0; lmax = mcget_ncount(); if (Flag_auto>0) Flag_auto=0; }
        if (!strcmp(token, "id") || !strcmp(token, "pixel"))
          { Set_Vars_Coord_Type = DEFS->COORD_PIXELID; 
            strcpy(Set_Vars_Coord_Label,"Pixel ID [1]"); 
            strcpy(Set_Vars_Coord_Var,"id"); lmin = 0; lmax = FLT_MAX; 
            if (Flag_auto>0) Flag_auto=0;
            Vars->Flag_List = 1; }
        if (!strcmp(token, "user") || !strcmp(token, "user1") || !strcmp(token, "u1"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER1; strncpy(Set_Vars_Coord_Label,Vars->UserName1,30); strcpy(Set_Vars_Coord_Var,"U1"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user2") || !strcmp(token, "u2"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER2; strncpy(Set_Vars_Coord_Label,Vars->UserName2,30); strcpy(Set_Vars_Coord_Var,"U2"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "user3") || !strcmp(token, "u3"))
          { Set_Vars_Coord_Type = DEFS->COORD_USER3; strncpy(Set_Vars_Coord_Label,Vars->UserName3,30); strcpy(Set_Vars_Coord_Var,"U3"); lmin = -1e10; lmax = 1e10; }

        if (!strcmp(token, "userdouble0") || !strcmp(token, "ud0"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE0; strcpy(Set_Vars_Coord_Label,"ud0 [1]"); strcpy(Set_Vars_Coord_Var,"ud0"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble1") || !strcmp(token, "ud1"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE1; strcpy(Set_Vars_Coord_Label,"ud1 [1]"); strcpy(Set_Vars_Coord_Var,"ud1"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble2") || !strcmp(token, "ud2"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE2; strcpy(Set_Vars_Coord_Label,"ud2 [1]"); strcpy(Set_Vars_Coord_Var,"ud2"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble3") || !strcmp(token, "ud3"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE3; strcpy(Set_Vars_Coord_Label,"ud3 [1]"); strcpy(Set_Vars_Coord_Var,"ud3"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble4") || !strcmp(token, "ud4"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE4; strcpy(Set_Vars_Coord_Label,"ud4 [1]"); strcpy(Set_Vars_Coord_Var,"ud4"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble5") || !strcmp(token, "ud5"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE5; strcpy(Set_Vars_Coord_Label,"ud5 [1]"); strcpy(Set_Vars_Coord_Var,"ud5"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble6") || !strcmp(token, "ud6"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE6; strcpy(Set_Vars_Coord_Label,"ud6 [1]"); strcpy(Set_Vars_Coord_Var,"ud6"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble7") || !strcmp(token, "ud7"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE7; strcpy(Set_Vars_Coord_Label,"ud7 [1]"); strcpy(Set_Vars_Coord_Var,"ud7"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble8") || !strcmp(token, "ud8"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE8; strcpy(Set_Vars_Coord_Label,"ud8 [1]"); strcpy(Set_Vars_Coord_Var,"ud8"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble9") || !strcmp(token, "ud9"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE9; strcpy(Set_Vars_Coord_Label,"ud9 [1]"); strcpy(Set_Vars_Coord_Var,"ud9"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble10") || !strcmp(token, "ud10"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE10; strcpy(Set_Vars_Coord_Label,"ud10 [1]"); strcpy(Set_Vars_Coord_Var,"ud10"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble11") || !strcmp(token, "ud11"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE11; strcpy(Set_Vars_Coord_Label,"ud11 [1]"); strcpy(Set_Vars_Coord_Var,"ud11"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble12") || !strcmp(token, "ud12"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE12; strcpy(Set_Vars_Coord_Label,"ud12 [1]"); strcpy(Set_Vars_Coord_Var,"ud12"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble13") || !strcmp(token, "ud13"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE13; strcpy(Set_Vars_Coord_Label,"ud13 [1]"); strcpy(Set_Vars_Coord_Var,"ud13"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble14") || !strcmp(token, "ud14"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE14; strcpy(Set_Vars_Coord_Label,"ud14 [1]"); strcpy(Set_Vars_Coord_Var,"ud14"); lmin = -1e10; lmax = 1e10; }
        if (!strcmp(token, "userdouble15") || !strcmp(token, "ud15"))
          { Set_Vars_Coord_Type = DEFS->COORD_USERDOUBLE15; strcpy(Set_Vars_Coord_Label,"ud15 [1]"); strcpy(Set_Vars_Coord_Var,"ud15"); lmin = -1e10; lmax = 1e10; }

        /* now stores variable keywords detected, if any */
        if (Set_Vars_Coord_Type != DEFS->COORD_NONE)
        {
          int Coord_Number = Vars->Coord_Number;
          if (Vars->Flag_log) { Set_Vars_Coord_Type |= DEFS->COORD_LOG; Vars->Flag_log = 0; }
          if (Flag_abs) { Set_Vars_Coord_Type |= DEFS->COORD_ABS; Flag_abs = 0; }
          if (Flag_auto != 0) { Set_Vars_Coord_Type |= DEFS->COORD_AUTO; 
            if (Flag_auto > 0) Flag_auto = 0; }
          if (Set_Coord_Mode == DEFS->COORD_SIGNAL)
          {
            Coord_Number = 0;
            Vars->Flag_signal = Set_Vars_Coord_Type;
          }
          else
          {
            if (Coord_Number < MONnD_COORD_NMAX)
            { Coord_Number++;
              Vars->Coord_Number = Coord_Number; 
              if (Set_Vars_Coord_Type != DEFS->COORD_PIXELID)
                Vars->Coord_NumberNoPixel++;
            }
            else if (Vars->Flag_Verbose) printf("Monitor_nD: %s reached max number of variables (%i).\n", Vars->compcurname, MONnD_COORD_NMAX);
          }
          Vars->Coord_Type[Coord_Number] = Set_Vars_Coord_Type;
          strncpy(Vars->Coord_Label[Coord_Number], Set_Vars_Coord_Label,30);
          strncpy(Vars->Coord_Var[Coord_Number], Set_Vars_Coord_Var,30);
          if (lmin > lmax) { XY = lmin; lmin=lmax; lmax = XY; }
          Vars->Coord_Min[Coord_Number] = lmin;
          Vars->Coord_Max[Coord_Number] = lmax;
          if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT || Set_Vars_Coord_Type == DEFS->COORD_PIXELID || Set_Vars_Coord_Type == DEFS->COORD_SIGNAL)
            Vars->Coord_Bin[Coord_Number] = 1;
          else
            Vars->Coord_Bin[Coord_Number] = 20;
          Set_Coord_Mode = DEFS->COORD_VAR;
          Flag_All = 0;
          Flag_No  = 0;
        } else {
          /* no variable name could be read from options */
          if (!iskeyword) {
            if (strcmp(token, "cm2") && strcmp(token, "incoming")
             && strcmp(token, "outgoing") && strcmp(token, "cm2")
             && strcmp(token, "cm^2") && strcmp(token, "float")
             && strcmp(token, "double") && strcmp(token, "binary")
             && strcmp(token, "steradian") && Vars->Flag_Verbose)
              printf("Monitor_nD: %s: unknown '%s' keyword in 'options'. Ignoring.\n", Vars->compcurname, token);
          }
        }
      carg++;
      } /* end if token */
    } /* end while carg */
    free(option_copy);
    if (carg == 128) printf("Monitor_nD: %s reached max number of tokens (%i). Skipping.\n", Vars->compcurname, 128);

    if ((Vars->Flag_Shape == DEFS->SHAPE_BOX) && (fabs(Vars->mzmax - Vars->mzmin) == 0)) Vars->Flag_Shape = DEFS->SHAPE_SQUARE;

    if (Vars->Flag_log == 1) Vars->Coord_Type[0] |= DEFS->COORD_LOG;
    if (Vars->Coord_Number == 0)
    { Vars->Flag_Auto_Limits=0; Vars->Flag_Multiple=0; Vars->Flag_List=0; }

    /* now setting Monitor Name from variable labels */
    strcpy(Vars->Monitor_Label,"");
    XY = 1; /* will contain total bin number */
    for (i = 0; i <= Vars->Coord_Number; i++)
    {
      if (Flag_auto != 0) Vars->Coord_Type[i] |= DEFS->COORD_AUTO;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if ((Set_Vars_Coord_Type == DEFS->COORD_X)
       || (Set_Vars_Coord_Type == DEFS->COORD_Y)
       || (Set_Vars_Coord_Type == DEFS->COORD_Z))
       strcpy(Short_Label[i],"Position");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_THETA)
       || (Set_Vars_Coord_Type == DEFS->COORD_PHI)
       || (Set_Vars_Coord_Type == DEFS->COORD_ANGLE))
       strcpy(Short_Label[i],"Angle");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_XY)
       || (Set_Vars_Coord_Type == DEFS->COORD_XZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_YZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_RADIUS))
       strcpy(Short_Label[i],"Radius");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_VX)
       || (Set_Vars_Coord_Type == DEFS->COORD_VY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_V)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_VYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_VXZ))
       strcpy(Short_Label[i],"Velocity");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_KX)
       || (Set_Vars_Coord_Type == DEFS->COORD_KY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXY)
       || (Set_Vars_Coord_Type == DEFS->COORD_KYZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_KXZ)
       || (Set_Vars_Coord_Type == DEFS->COORD_K))
       strcpy(Short_Label[i],"Wavevector");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_SX)
       || (Set_Vars_Coord_Type == DEFS->COORD_SY)
       || (Set_Vars_Coord_Type == DEFS->COORD_SZ))
       strcpy(Short_Label[i],"Spin");
      else
      if ((Set_Vars_Coord_Type == DEFS->COORD_HDIV)
       || (Set_Vars_Coord_Type == DEFS->COORD_VDIV))
       strcpy(Short_Label[i],"Divergence");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY)
       strcpy(Short_Label[i],"Energy");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA)
       strcpy(Short_Label[i],"Wavelength");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT)
       strcpy(Short_Label[i],"Neutron_ID");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID)
       strcpy(Short_Label[i],"Pixel_ID");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_T)
          strcpy(Short_Label[i],"Time_Of_Flight");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_P)
          strcpy(Short_Label[i],"Intensity");
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER1)
          strncpy(Short_Label[i],Vars->UserName1,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER2)
          strncpy(Short_Label[i],Vars->UserName2,30);
      else
      if (Set_Vars_Coord_Type == DEFS->COORD_USER3)
          strncpy(Short_Label[i],Vars->UserName3,30);
      else
          strcpy(Short_Label[i],"Unknown");

      if (Vars->Coord_Type[i] & DEFS->COORD_ABS)
      { strcat(Vars->Coord_Label[i]," (abs)"); }

      if (Vars->Coord_Type[i] & DEFS->COORD_LOG)
      { strcat(Vars->Coord_Label[i]," (log)"); }

      strcat(Vars->Monitor_Label, " ");
      strcat(Vars->Monitor_Label, Short_Label[i]);
      XY *= Vars->Coord_Bin[i];

    } /* end for Short_Label */

    if ((Vars->Coord_Type[0] & (DEFS->COORD_LOG-1)) == DEFS->COORD_P) {
      strncat(Vars->Coord_Label[0], " [n/s", 30);
      if (Vars->Flag_per_cm2) strncat(Vars->Coord_Label[0], "/cm2", 30);

      if (XY > 1 && Vars->Coord_Number)
        strncat(Vars->Coord_Label[0], "/bin", 30);
      strncat(Vars->Coord_Label[0], "]", 30);
    }

    /* update label 'signal per bin' if more than 1 bin */
    if (XY > 1 && Vars->Coord_Number) {
      if (Vars->Flag_capture)
        printf("Monitor_nD: %s: Using capture flux weightening on %ld bins.\n"
               "WARNING     Use binned data with caution, and prefer monitor integral value (I,Ierr).\n", Vars->compcurname, (long)XY);
    }

    strcat(Vars->Monitor_Label, " Monitor");
    if (Vars->Flag_Shape == DEFS->SHAPE_SQUARE) strcat(Vars->Monitor_Label, " (Square)");
    if (Vars->Flag_Shape == DEFS->SHAPE_DISK)   strcat(Vars->Monitor_Label, " (Disk)");
    if (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) strcat(Vars->Monitor_Label, " (Sphere)");
    if (Vars->Flag_Shape == DEFS->SHAPE_CYLIND) strcat(Vars->Monitor_Label, " (Cylinder)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BANANA) strcat(Vars->Monitor_Label, " (Banana)");
    if (Vars->Flag_Shape == DEFS->SHAPE_BOX)    strcat(Vars->Monitor_Label, " (Box)");
    if (Vars->Flag_Shape == DEFS->SHAPE_PREVIOUS) strcat(Vars->Monitor_Label, " (on PREVIOUS)");
    if (Vars->Flag_Shape == DEFS->SHAPE_OFF) strcat(Vars->Monitor_Label, " (OFF geometry)");
    if ((Vars->Flag_Shape == DEFS->SHAPE_CYLIND) || (Vars->Flag_Shape == DEFS->SHAPE_BANANA) || (Vars->Flag_Shape == DEFS->SHAPE_SPHERE) || (Vars->Flag_Shape == DEFS->SHAPE_BOX))
    {
      if (strstr(Vars->option, "incoming"))
      {
        Vars->Flag_Shape = abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [in]");
      }
      else /* if strstr(Vars->option, "outgoing")) */
      {
        Vars->Flag_Shape = -abs(Vars->Flag_Shape);
        strcat(Vars->Monitor_Label, " [out]");
      }
    }
    if (Vars->Flag_UsePreMonitor == 1)
    {
        strcat(Vars->Monitor_Label, " at ");
        strncat(Vars->Monitor_Label, Vars->UserName1,30);
    }
    if (Vars->Flag_log == 1) strcat(Vars->Monitor_Label, " [log] ");

    /* now allocate memory to store variables in TRACE */

    /* Vars->Coord_Number  0   : intensity or signal
     * Vars->Coord_Number  1:n : detector variables */

    if ((Vars->Coord_NumberNoPixel != 2) && !Vars->Flag_Multiple && !Vars->Flag_List)
    { Vars->Flag_Multiple = 1; /* default is n1D */
      if (Vars->Coord_Number != Vars->Coord_NumberNoPixel) Vars->Flag_List = 1; }

    /* list and auto limits case : Vars->Flag_List or Vars->Flag_Auto_Limits
     * -> Buffer to flush and suppress after Vars->Flag_Auto_Limits
     */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      Vars->Mon2D_Buffer = (double *)malloc((Vars->Coord_Number+1)*Vars->Buffer_Block*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
      { printf("Monitor_nD: %s cannot allocate Vars->Mon2D_Buffer (%li). No list and auto limits.\n", Vars->compcurname, Vars->Buffer_Block*(Vars->Coord_Number+1)*sizeof(double)); Vars->Flag_List = 0; Vars->Flag_Auto_Limits = 0; }
      else
      {
        for (i=0; i < (Vars->Coord_Number+1)*Vars->Buffer_Block; Vars->Mon2D_Buffer[i++] = (double)0);
      }
      Vars->Buffer_Size = Vars->Buffer_Block;
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_NumberNoPixel)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Number)*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, (Vars->Coord_Number)*sizeof(double *)); exit(-1); }
      for (i= 1; i <= Vars->Coord_Number; i++)
      {
        Vars->Mon2D_N[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p[i-1]  = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        Vars->Mon2D_p2[i-1] = (double *)malloc(Vars->Coord_Bin[i]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s n1D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[i], i, (Vars->Coord_Bin[i])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[i]; j++ )
          { Vars->Mon2D_N[i-1][j] = (double)0; Vars->Mon2D_p[i-1][j] = (double)0; Vars->Mon2D_p2[i-1][j] = (double)0; }
        }
      }
    }
    else /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_NumberNoPixel == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      Vars->Mon2D_N  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p  = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      Vars->Mon2D_p2 = (double **)malloc((Vars->Coord_Bin[1])*sizeof(double *));
      if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
      { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2 (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], (Vars->Coord_Bin[1])*sizeof(double *)); exit(-1); }
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        Vars->Mon2D_N[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p[i]  = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        Vars->Mon2D_p2[i] = (double *)malloc(Vars->Coord_Bin[2]*sizeof(double));
        if ((Vars->Mon2D_N == NULL) || (Vars->Mon2D_p == NULL) || (Vars->Mon2D_p2 == NULL))
        { fprintf(stderr,"Monitor_nD: %s 2D cannot allocate %s Vars->Mon2D_N/p/p2[%li] (%li). Fatal.\n", Vars->compcurname, Vars->Coord_Var[1], i, (Vars->Coord_Bin[2])*sizeof(double *)); exit(-1); }
        else
        {
          for (j=0; j < Vars->Coord_Bin[2]; j++ )
          { Vars->Mon2D_N[i][j] = (double)0; Vars->Mon2D_p[i][j] = (double)0; Vars->Mon2D_p2[i][j] = (double)0; }
        }
      }
    }
    else {
      Vars->Mon2D_N = Vars->Mon2D_p = Vars->Mon2D_p2 = NULL;
    }
      /* no Mon2D allocated for
       * (Vars->Coord_Number != 2) && !Vars->Flag_Multiple && Vars->Flag_List */

    Vars->psum  = 0;
    Vars->p2sum = 0;
    Vars->Nsum  = 0;

    Vars->area  = fabs(Vars->mxmax - Vars->mxmin)*fabs(Vars->mymax - Vars->mymin)*1E4; /* in cm**2 for square and box shapes */
    Vars->Sphere_Radius = fabs(Vars->mxmax - Vars->mxmin)/2;
    if ((abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    {
      Vars->area = PI*Vars->Sphere_Radius*Vars->Sphere_Radius*1E4; /* disk shapes */
    }


    if (Vars->area == 0 && abs(Vars->Flag_Shape) != DEFS->SHAPE_PREVIOUS ) {
      if (abs(Vars->Flag_Shape) != DEFS->SHAPE_OFF) {  
	Vars->Coord_Number = 0;
      }
    }
    if (Vars->Coord_Number == 0 && Vars->Flag_Verbose)
      printf("Monitor_nD: %s is unactivated (0D)\n", Vars->compcurname);
    Vars->Cylinder_Height = fabs(Vars->mymax - Vars->mymin);

    if (Vars->Flag_Verbose)
    {
      printf("Monitor_nD: %s is a %s.\n", Vars->compcurname, Vars->Monitor_Label);
      printf("Monitor_nD: version %s with options=%s\n", MONITOR_ND_LIB_H, Vars->option);
    }
    
    /* compute the product of bin dimensions for PixelID */
    Vars->Coord_BinProd[0]=1;
    for (i = 1; i <= Vars->Coord_Number; i++)
      Vars->Coord_BinProd[i]=Vars->Coord_Bin[i]*Vars->Coord_BinProd[i-1];
  } /* end Monitor_nD_Init */

/* ========================================================================= */
/* Monitor_nD_Trace: this routine is used to monitor one propagating neutron */
/* return values: 0=neutron was absorbed, -1=neutron was outside bounds, 1=neutron was measured*/
/* ========================================================================= */

int Monitor_nD_Trace(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars, _class_particle* _particle)
{

  double  XY=0, pp=0;
  long    i =0, j =0;
  double  Coord[MONnD_COORD_NMAX];
  long    Coord_Index[MONnD_COORD_NMAX];
  char    While_End   =0;
  long    While_Buffer=0;
  char    Set_Vars_Coord_Type = DEFS->COORD_NONE;
  
  /* the logic below depends mainly on:
       Flag_List:        1=store 1 buffer, 2=list all, 3=re-use buffer 
       Flag_Auto_Limits: 0 (no auto limits/list), 1 (store events into Buffer), 2 (re-emit store events)
   */

  /* Vars->Flag_Auto_Limits=1: buffer full, we read the Buffer, and determine min and max bounds */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 1) && (Vars->Coord_Number > 0))
  {
    /* auto limits case : get limits in Buffer for each variable */
          /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
    if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li events) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
    for (i = 1; i <= Vars->Coord_Number; i++)
    {
      if (Vars->Coord_Type[i] & DEFS->COORD_AUTO)
      {
        Vars->Coord_Min[i] =  FLT_MAX;
        Vars->Coord_Max[i] = -FLT_MAX;
        for (j = 0; j < Vars->Buffer_Counter; j++)
        {
          XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
          if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
          if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
        }
        if  (Vars->Flag_Verbose)  
          printf("  %s: min=%g max=%g\n", Vars->Coord_Var[i], Vars->Coord_Min[i], Vars->Coord_Max[i]);
      }
    }
    Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step (read Buffer and generate new events to store in histograms) */
  } /* end if Flag_Auto_Limits == 1 */

#ifndef OPENACC
  /* manage realloc for 'list all' if Buffer size exceeded: flush Buffer to file */
  if ((Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List >= 2))
  {
    if (Vars->Buffer_Size >= 1000000 || Vars->Flag_List == 3)
    { /* save current (possibly append) and re-use Buffer */

      Monitor_nD_Save(DEFS, Vars);
      Vars->Flag_List = 3;
      Vars->Buffer_Block = Vars->Buffer_Size;
      Vars->Buffer_Counter  = 0;
      Vars->Neutron_Counter = 0;

    }
    else
    {
      Vars->Mon2D_Buffer  = (double *)realloc(Vars->Mon2D_Buffer, (Vars->Coord_Number+1)*(Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double));
      if (Vars->Mon2D_Buffer == NULL)
            { printf("Monitor_nD: %s cannot reallocate Vars->Mon2D_Buffer[%li] (%li). Skipping.\n", Vars->compcurname, i, (long int)(Vars->Neutron_Counter+Vars->Buffer_Block)*sizeof(double)); Vars->Flag_List = 1; }
      else { Vars->Buffer_Counter = 0; Vars->Buffer_Size = Vars->Neutron_Counter+Vars->Buffer_Block; }
    }
  } /* end if Buffer realloc */
#endif

  char    outsidebounds=0;
  while (!While_End)
  { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) or passing neutron */
    if ((Vars->Flag_Auto_Limits == 2) && (Vars->Coord_Number > 0))
    { /* Vars->Flag_Auto_Limits == 2: read back from Buffer (Buffer is filled or auto limits have been computed) */
      if (While_Buffer < Vars->Buffer_Block)
      {
        /* first while loop (While_Buffer) */
        /* auto limits case : scan Buffer within limits and store in Mon2D */
        Coord[0] = pp = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

        for (i = 1; i <= Vars->Coord_Number; i++)
        {
          /* scanning variables in Buffer */
          if (Vars->Coord_Bin[i] <= 1) continue;
          XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);

          Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
          if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
          else        Coord_Index[i] = 0;
          if (Vars->Flag_With_Borders)
          {
            if (Coord_Index[i] < 0)                   Coord_Index[i] = 0;
            if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
          }
        } /* end for */
        
        /* update the PixelID, we compute it from the previous variables index */
        if (Vars->Coord_NumberNoPixel < Vars->Coord_Number) /* there is a Pixel variable */
        for (i = 1; i <= Vars->Coord_Number; i++) {
          char Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
          if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID) {
            char flag_outside=0;
            Coord_Index[i] = Coord[i] = 0;
            for (j= 1; j < i; j++) {
              /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
              if (Vars->Coord_Bin[j] == 1) continue; 
              if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) {
                flag_outside=1;
                Coord[i] = 0;
                break;
              }
              Coord[i] += Coord_Index[j]*Vars->Coord_BinProd[j-1];
            }
            if (!flag_outside) {
              Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)] = Coord[i];
            }
          } /* end if PixelID */
        }
        While_Buffer++;
      } /* end if in Buffer */
      else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
      {
        Vars->Flag_Auto_Limits = 0;
        if (!Vars->Flag_List) /* free Buffer not needed anymore (no list to output) */
        { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, p2) */
          free(Vars->Mon2D_Buffer); Vars->Mon2D_Buffer = NULL;
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li) in TRACE.\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
      }
    } /* if Vars->Flag_Auto_Limits == 2 */
    
    if (Vars->Flag_Auto_Limits != 2 || !Vars->Coord_Number) /* Vars->Flag_Auto_Limits == 0 (no auto limits/list) or 1 (store events into Buffer) */
    {
      /* automatically compute area and steradian solid angle when in AUTO mode */
      /* compute the steradian solid angle incoming on the monitor */
      double v;
      double tmp;
      v=sqrt(_particle->vx*_particle->vx + _particle->vy*_particle->vy + _particle->vz*_particle->vz);
      tmp=_particle->x;
      if (Vars->min_x > _particle->x){
        #pragma acc atomic write
        Vars->min_x = tmp;
      }
      if (Vars->max_x < _particle->x){
        #pragma acc atomic write
        Vars->max_x = tmp;
      }
      tmp=_particle->y;
      if (Vars->min_y > _particle->y){
        #pragma acc atomic write
        Vars->min_y = tmp;
      }
      if (Vars->max_y < _particle->y){
	tmp=_particle->y;
        #pragma acc atomic write
	Vars->max_y = tmp;
      }

      #pragma acc atomic
      Vars->mean_p = Vars->mean_p + _particle->p;
      if (v) {
        tmp=_particle->p*fabs(_particle->vx/v);
        #pragma acc atomic
        Vars->mean_dx = Vars->mean_dx + tmp; //_particle->p*fabs(_particle->vx/v);
        tmp=_particle->p*fabs(_particle->vy/v);
        #pragma acc atomic
        Vars->mean_dy = Vars->mean_dy + tmp; //_particle->p*fabs(_particle->vy/v);
      }

      for (i = 0; i <= Vars->Coord_Number; i++)
      { /* handle current neutron : last while */
        XY = 0;
        Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
        /* get values for variables to monitor */
        if (Set_Vars_Coord_Type == DEFS->COORD_X) XY = _particle->x;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Y) XY = _particle->y;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_Z) XY = _particle->z;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VX) XY = _particle->vx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VY) XY = _particle->vy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VZ) XY = _particle->vz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KX) XY = V2K*_particle->vx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KY) XY = V2K*_particle->vy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KZ) XY = V2K*_particle->vz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SX) XY = _particle->sx;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SY) XY = _particle->sy;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_SZ) XY = _particle->sz;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_T) XY = _particle->t;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_P) XY = _particle->p;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE0) XY = Vars->UserDoubles[0];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE1) XY = Vars->UserDoubles[1];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE2) XY = Vars->UserDoubles[2];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE3) XY = Vars->UserDoubles[3];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE4) XY = Vars->UserDoubles[4];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE5) XY = Vars->UserDoubles[5];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE6) XY = Vars->UserDoubles[6];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE7) XY = Vars->UserDoubles[7];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE8) XY = Vars->UserDoubles[8];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE9) XY = Vars->UserDoubles[9];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE10) XY = Vars->UserDoubles[10];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE11) XY = Vars->UserDoubles[11];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE12) XY = Vars->UserDoubles[12];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE13) XY = Vars->UserDoubles[13];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE14) XY = Vars->UserDoubles[14];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USERDOUBLE15) XY = Vars->UserDoubles[15];
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_HDIV) XY = RAD2DEG*atan2(_particle->vx,_particle->vz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VDIV) XY = RAD2DEG*atan2(_particle->vy,_particle->vz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_V) XY = sqrt(_particle->vx*_particle->vx+_particle->vy*_particle->vy+_particle->vz*_particle->vz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
          XY = sqrt(_particle->x*_particle->x+_particle->y*_particle->y+_particle->z*_particle->z);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XY)
          XY = sqrt(_particle->x*_particle->x+_particle->y*_particle->y)*(_particle->x > 0 ? 1 : -1);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_YZ) XY = sqrt(_particle->y*_particle->y+_particle->z*_particle->z);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_XZ)
          XY = sqrt(_particle->x*_particle->x+_particle->z*_particle->z);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXY) XY = sqrt(_particle->vx*_particle->vx+_particle->vy*_particle->vy);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VXZ) XY = sqrt(_particle->vx*_particle->vx+_particle->vz*_particle->vz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_VYZ) XY = sqrt(_particle->vy*_particle->vy+_particle->vz*_particle->vz);
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_K) { XY = sqrt(_particle->vx*_particle->vx+_particle->vy*_particle->vy+_particle->vz*_particle->vz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXY) { XY = sqrt(_particle->vx*_particle->vx+_particle->vy*_particle->vy);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KXZ) { XY = sqrt(_particle->vx*_particle->vx+_particle->vz*_particle->vz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_KYZ) { XY = sqrt(_particle->vy*_particle->vy+_particle->vz*_particle->vz);  XY *= V2K; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ENERGY) { XY = _particle->vx*_particle->vx+_particle->vy*_particle->vy+_particle->vz*_particle->vz;  XY *= VS2E; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_LAMBDA) { XY = sqrt(_particle->vx*_particle->vx+_particle->vy*_particle->vy+_particle->vz*_particle->vz);  XY *= V2K; if (XY != 0) XY = 2*PI/XY; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_NCOUNT) XY = Vars->Neutron_Counter;
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
        {  XY = sqrt(_particle->vx*_particle->vx+_particle->vy*_particle->vy);
           if (_particle->vz != 0)
                XY = RAD2DEG*atan2(XY,_particle->vz)*(_particle->x > 0 ? 1 : -1);
           else XY = 0;
        }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_THETA)  { if (_particle->z != 0) XY = RAD2DEG*atan2(_particle->x,_particle->z); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PHI) { double rr=sqrt(_particle->x*_particle->x+ _particle->y*_particle->y + _particle->z*_particle->z); if (rr != 0) XY = RAD2DEG*asin(_particle->y/rr); }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER1) {int fail; XY = particle_getvar(_particle,Vars->UserVariable1,&fail); if(fail) XY=0; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER2) {int fail; XY = particle_getvar(_particle,Vars->UserVariable2,&fail); if(fail) XY=0; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_USER3) {int fail; XY = particle_getvar(_particle,Vars->UserVariable3,&fail); if(fail) XY=0; }
        else
        if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID && !Vars->Flag_Auto_Limits) {
          /* compute the PixelID from previous coordinates 
             the PixelID is the product of Coord_Index[i] in the detector geometry 
             pixelID = sum( Coord_Index[j]*prod(Vars->Coord_Bin[1:(j-1)]) )
             
             this does not apply when we store events in the buffer as Coord_Index
             is not set. Then the pixelID will be re-computed during SAVE.
          */
          char flag_outside=0;
          for (j= 1; j < i; j++) {
            /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
            if (Vars->Coord_Bin[j] <= 1) continue; 
            if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) { 
              flag_outside=1; XY=0; break;
            }
            XY += Coord_Index[j]*Vars->Coord_BinProd[j-1];
          }
	  if (Vars->Flag_mantid && Vars->Flag_OFF && Vars->OFF_polyidx >=0) XY=Vars->OFF_polyidx;
          if (!flag_outside) XY += Vars->Coord_Min[i];
        }
        
        /* handle 'abs' and 'log' keywords */
        if (Vars->Coord_Type[i] & DEFS->COORD_ABS) XY=fabs(XY);

        if (Vars->Coord_Type[i] & DEFS->COORD_LOG) /* compute log of variable if requested */
        {  if (XY > 0) XY = log(XY)/log(10);
           else        XY = -100; }

        Coord[i] = XY; Coord_Index[i] = 0;
        if (i == 0) { pp = XY; Coord_Index[i] = 0; }
        else {
        /* check bounds for variables which have no automatic limits */
          if ((!Vars->Flag_Auto_Limits || !(Vars->Coord_Type[i] & DEFS->COORD_AUTO)) && Vars->Coord_Bin[i]>1)
          { /* compute index in histograms for each variable to monitor */
            XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
            if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
            if (Vars->Flag_With_Borders)
            {
              if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
              if (Coord_Index[i] < 0) Coord_Index[i] = 0;
            }
            //if (0 > Coord_Index[i] || Coord_Index[i] >= Vars->Coord_Bin[i])
            //  outsidebounds=1;
          } /* else will get Index later from Buffer when Flag_Auto_Limits == 2 */
        }
        
      } /* end for i */
      While_End = 1;
    }/* end else if Vars->Flag_Auto_Limits == 2 */
    
    /* ====================================================================== */
    /* store n1d/2d neutron from Buffer (Auto_Limits == 2) or current neutron in while */
    if (Vars->Flag_Auto_Limits != 1) /* not when storing auto limits Buffer */
    {
      /* apply per cm2 */
      if (Vars->Flag_per_cm2 && Vars->area != 0)
        pp /= Vars->area;

      /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
      if ( Vars->Coord_NumberNoPixel == 2 && !Vars->Flag_Multiple)
      { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
        
        i = Coord_Index[1];
        j = Coord_Index[2];
        if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
        {
          if (Vars->Mon2D_N) {
	    double p2 = pp*pp;
            #pragma acc atomic
	    Vars->Mon2D_N[i][j] = Vars->Mon2D_N[i][j]+1;
            #pragma acc atomic
	    Vars->Mon2D_p[i][j] = Vars->Mon2D_p[i][j]+pp;
            #pragma acc atomic
	    Vars->Mon2D_p2[i][j] = Vars->Mon2D_p2[i][j] + p2;
	  }
        } else {
          outsidebounds=1; 
        }
      } else {
        /* 1D and n1D case : Vars->Flag_Multiple */
        /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          
        for (i= 1; i <= Vars->Coord_Number; i++) {
          j = Coord_Index[i];
          if (j >= 0 && j < Vars->Coord_Bin[i]) {
            if  (Vars->Flag_Multiple && Vars->Mon2D_N) {
	      if (Vars->Mon2D_N) {
		double p2 = pp*pp;
                #pragma acc atomic
		Vars->Mon2D_N[i-1][j] = Vars->Mon2D_N[i-1][j]+1;
                #pragma acc atomic
		Vars->Mon2D_p[i-1][j] = Vars->Mon2D_p[i-1][j]+pp;
		#pragma acc atomic
		Vars->Mon2D_p2[i-1][j] = Vars->Mon2D_p2[i-1][j] + p2;
	      }
	    }
          } else { 
            outsidebounds=1;
            break;
          }
        }
      }
    } /* end (Vars->Flag_Auto_Limits != 1) */
    
    if (Vars->Flag_Auto_Limits != 2 && !outsidebounds) /* not when reading auto limits Buffer */
    { /* now store Coord into Buffer (no index needed) if necessary (list or auto limits) */
      if ((Vars->Buffer_Counter < Vars->Buffer_Block) && ((Vars->Flag_List) || (Vars->Flag_Auto_Limits == 1)))
      {
        #pragma acc atomic
        Vars->Neutron_Counter = Vars->Neutron_Counter + 1;
        for (i = 0; i <= Vars->Coord_Number; i++)
        {
	  // This is is where the list is appended. How to make this "atomic"?
          #pragma acc atomic write 
          Vars->Mon2D_Buffer[i + Vars->Neutron_Counter*(Vars->Coord_Number+1)] = Coord[i];
        }
	#pragma acc atomic 
        Vars->Buffer_Counter = Vars->Buffer_Counter + 1;
        if (Vars->Flag_Verbose && (Vars->Buffer_Counter >= Vars->Buffer_Block) && (Vars->Flag_List == 1)) 
          printf("Monitor_nD: %s %li neutrons stored in List.\n", Vars->compcurname, Vars->Buffer_Counter);
      }
    } /* end (Vars->Flag_Auto_Limits != 2) */
    
  } /* end while */
  #pragma acc atomic
  Vars->Nsum = Vars->Nsum + 1;
  #pragma acc atomic
  Vars->psum  = Vars->psum + pp;
  #pragma acc atomic
  Vars->p2sum = Vars->p2sum + pp*pp;

  /*determine return value: 1:neutron was in bounds and measured, -1: outside bounds, 0: outside bounds, should be absorbed.*/
  if(outsidebounds){
      if(Vars->Flag_Absorb){
          return 0;
      }else{
          return -1;
      }
  }
  return 1;
} /* end Monitor_nD_Trace */

/* ========================================================================= */
/* Monitor_nD_Save: this routine is used to save data files                  */
/* ========================================================================= */

MCDETECTOR Monitor_nD_Save(MonitornD_Defines_type *DEFS, MonitornD_Variables_type *Vars)
  {
    char   *fname;
    long    i,j;
    double *p0m = NULL;
    double *p1m = NULL;
    double *p2m = NULL;
    char    Coord_X_Label[CHAR_BUF_LENGTH];
    double  min1d, max1d;
    double  min2d, max2d;
    char    While_End = 0;
    long    While_Buffer = 0;
    double  XY=0, pp=0;
    double  Coord[MONnD_COORD_NMAX];
    long    Coord_Index[MONnD_COORD_NMAX];
    char    label[CHAR_BUF_LENGTH];

    MCDETECTOR detector;

    if (Vars->Flag_Verbose && Vars->Flag_per_cm2) {
      printf("Monitor_nD: %s: active flat detector area is %g [cm^2], total area is %g [cm^2]\n",
        Vars->compcurname, (Vars->max_x-Vars->min_x)
                          *(Vars->max_y-Vars->min_y)*1E4, Vars->area);
      printf("Monitor_nD: %s: beam solid angle is %g [st] (%g x %g [deg^2])\n",
        Vars->compcurname,
        2*fabs(2*atan(Vars->mean_dx/Vars->mean_p)
         *sin(2*atan(Vars->mean_dy/Vars->mean_p)/2)),
        atan(Vars->mean_dx/Vars->mean_p)*RAD2DEG,
        atan(Vars->mean_dy/Vars->mean_p)*RAD2DEG);
    }

    /* check Buffer flush when end of simulation reached */
    if ((Vars->Buffer_Counter <= Vars->Buffer_Block) && Vars->Flag_Auto_Limits && Vars->Mon2D_Buffer && Vars->Buffer_Counter)
    {
      /* Get Auto Limits */
      if (Vars->Flag_Verbose) printf("Monitor_nD: %s getting %li Auto Limits from List (%li events).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);

      for (i = 1; i <= Vars->Coord_Number; i++)
      {
        if ((Vars->Coord_Type[i] & DEFS->COORD_AUTO) && Vars->Coord_Bin[i] > 1)
        {
          Vars->Coord_Min[i] = FLT_MAX;
          Vars->Coord_Max[i] = -FLT_MAX;
          for (j = 0; j < Vars->Buffer_Counter; j++)
          {
            XY = Vars->Mon2D_Buffer[i+j*(Vars->Coord_Number+1)];  /* scanning variables in Buffer */
            if (XY < Vars->Coord_Min[i]) Vars->Coord_Min[i] = XY;
            if (XY > Vars->Coord_Max[i]) Vars->Coord_Max[i] = XY;
          }
          if  (Vars->Flag_Verbose)  
            printf("  %s: min=%g max=%g in %li bins\n", Vars->Coord_Var[i], Vars->Coord_Min[i], Vars->Coord_Max[i], Vars->Coord_Bin[i]);
        }
      }
      Vars->Flag_Auto_Limits = 2;  /* pass to 2nd auto limits step */
      Vars->Buffer_Block = Vars->Buffer_Counter;

      while (!While_End)
      { /* we generate Coord[] and Coord_index[] from Buffer (auto limits) */
        /* simulation ended before Buffer was filled. Limits have to be computed, and stored events must be sent into histograms */
        
        if (While_Buffer < Vars->Buffer_Block)
        {
          /* first while loops (While_Buffer) */
          Coord[0] = Vars->Mon2D_Buffer[While_Buffer*(Vars->Coord_Number+1)];

          /* auto limits case : scan Buffer within limits and store in Mon2D */
          for (i = 1; i <= Vars->Coord_Number; i++)
          {
            /* scanning variables in Buffer */
            if (Vars->Coord_Bin[i] <= 1) Coord_Index[i] = 0;
            else {
              XY = (Vars->Coord_Max[i]-Vars->Coord_Min[i]);
              Coord[i] = Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)];
              if (XY > 0) Coord_Index[i] = floor((Coord[i]-Vars->Coord_Min[i])*Vars->Coord_Bin[i]/XY);
              else Coord_Index[i] = 0;
              if (Vars->Flag_With_Borders)
              {
                if (Coord_Index[i] < 0) Coord_Index[i] = 0;
                if (Coord_Index[i] >= Vars->Coord_Bin[i]) Coord_Index[i] = Vars->Coord_Bin[i] - 1;
              }
            }
          } /* end for */

          /* update the PixelID, we compute it from the previous variables index */
          for (i = 1; i <= Vars->Coord_Number; i++) {
            char Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
            if (Set_Vars_Coord_Type == DEFS->COORD_PIXELID) {
              char outsidebounds=0;
              Coord_Index[i] = Coord[i] = 0;
              for (j= 1; j < i; j++) {
                /* not for 1D variables with Bin=1 such as PixelID, NCOUNT, Intensity */
                if (Vars->Coord_Bin[j] == 1) continue; 
                if (0 > Coord_Index[j] || Coord_Index[j] >= Vars->Coord_Bin[j]) {
                  outsidebounds=1;
                  Coord[i] = 0;
                  break;
                }
                Coord[i] += Coord_Index[j]*Vars->Coord_BinProd[j-1];
              }
              if (!outsidebounds) {
                Vars->Mon2D_Buffer[i+While_Buffer*(Vars->Coord_Number+1)] = Coord[i];
              }
            } /* end if PixelID */
          }
          While_Buffer++;
        } /* end if in Buffer */
        else /* (While_Buffer >= Vars->Buffer_Block) && (Vars->Flag_Auto_Limits == 2) */
        {
          Vars->Flag_Auto_Limits = 0;
          While_End = 1;
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s flushed %li Auto Limits from List (%li).\n", Vars->compcurname, Vars->Coord_Number, Vars->Buffer_Counter);
        }

        /* store n1d/2d section from Buffer */

        pp = Coord[0];
        /* apply per cm2 or per st */
        if (Vars->Flag_per_cm2 && Vars->area      != 0)
          pp /= Vars->area;
        
        /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
        if (!Vars->Flag_Multiple && Vars->Coord_NumberNoPixel == 2)
        { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
          i = Coord_Index[1];
          j = Coord_Index[2];
          if (i >= 0 && i < Vars->Coord_Bin[1] && j >= 0 && j < Vars->Coord_Bin[2])
          {
            if (Vars->Mon2D_N) {
              Vars->Mon2D_N[i][j]++;
              Vars->Mon2D_p[i][j] += pp;
              Vars->Mon2D_p2[i][j] += pp*pp;
            }
          } else if (Vars->Flag_Absorb) pp=0;
        }
        else
        /* 1D and n1D case : Vars->Flag_Multiple */
        { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors (intensity is not included) */
          for (i= 1; i <= Vars->Coord_Number; i++)
          {
            j = Coord_Index[i];
            if (j >= 0 && j < Vars->Coord_Bin[i])
            {
              if (Vars->Flag_Multiple && Vars->Mon2D_N) {
                Vars->Mon2D_N[i-1][j]++;
                Vars->Mon2D_p[i-1][j] += pp;
                Vars->Mon2D_p2[i-1][j] += pp*pp;
              }
            } else if (Vars->Flag_Absorb) {
              pp=0; break;
            }
          }
        } /* end store 2D/1D */
        
      } /* end while */
    } /* end Force Get Limits */

    /* write output files (sent to file as p[i*n + j] vectors) */
    if (Vars->Coord_Number == 0)
    {
      double Nsum;
      double psum, p2sum;
      Nsum = Vars->Nsum;
      psum = Vars->psum;
      p2sum= Vars->p2sum;
      if (Vars->Flag_signal != DEFS->COORD_P && Nsum > 0)
      { psum /=Nsum; p2sum /= Nsum*Nsum; }
      /* DETECTOR_OUT_0D(Vars->Monitor_Label, Vars->Nsum, Vars->psum, Vars->p2sum); */
      detector = mcdetector_out_0D(Vars->Monitor_Label, Nsum, psum, p2sum, Vars->compcurname, Vars->compcurpos);
    }
    else
    if (strlen(Vars->Mon_File) > 0)
    {
      fname = (char*)malloc(strlen(Vars->Mon_File)+10*Vars->Coord_Number);
      if (Vars->Flag_List && Vars->Mon2D_Buffer) /* List: DETECTOR_OUT_2D */
      {
       
        if (Vars->Flag_List >= 2) Vars->Buffer_Size = Vars->Neutron_Counter;
        if (Vars->Buffer_Size >= Vars->Neutron_Counter)
          Vars->Buffer_Size = Vars->Neutron_Counter;
        strcpy(fname,Vars->Mon_File);
        if (strchr(Vars->Mon_File,'.') == NULL) strcat(fname, "_list");

        strcpy(Coord_X_Label,"");
        for (i= 0; i <= Vars->Coord_Number; i++)
        {
          strcat(Coord_X_Label, Vars->Coord_Var[i]);
          strcat(Coord_X_Label, " ");
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i]); }
        }
        if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s List (%lix%li).\n", Vars->compcurname, fname,(long int)Vars->Neutron_Counter,Vars->Coord_Number);

        /* handle the type of list output */
        strcpy(label, Vars->Monitor_Label);
        
        detector = mcdetector_out_list(
              label, "List of neutron events", Coord_X_Label,
              -Vars->Buffer_Size, Vars->Coord_Number+1,
              Vars->Mon2D_Buffer,
              fname, Vars->compcurname, Vars->compcurpos);
      }
      if (Vars->Flag_Multiple) /* n1D: DETECTOR_OUT_1D */
      {
        for (i= 0; i < Vars->Coord_Number; i++)
        {

          strcpy(fname,Vars->Mon_File);
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[i+1]); }
          sprintf(Coord_X_Label, "%s monitor", Vars->Coord_Label[i+1]);
          strcpy(label, Coord_X_Label);
          if (Vars->Coord_Bin[i+1] > 0) { /* 1D monitor */
            if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 1D (%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[i+1]);
            min1d = Vars->Coord_Min[i+1];
            max1d = Vars->Coord_Max[i+1];
            if (min1d == max1d) max1d = min1d+1e-6;
            p1m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            p2m = (double *)malloc(Vars->Coord_Bin[i+1]*sizeof(double));
            if (p2m == NULL) /* use Raw Buffer line output */
            {
              if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for output. Using raw data.\n", Vars->compcurname);
              if (p1m != NULL) free(p1m);
              detector = mcdetector_out_1D(
              label,
              Vars->Coord_Label[i+1],
              Vars->Coord_Label[0],
              Vars->Coord_Var[i+1],
              min1d, max1d,
              Vars->Coord_Bin[i+1],
              Vars->Mon2D_N[i],Vars->Mon2D_p[i],Vars->Mon2D_p2[i],
              fname, Vars->compcurname, Vars->compcurpos);
            } /* if (p2m == NULL) */
            else
            {
              if (Vars->Flag_log != 0)
              {
                XY = FLT_MAX;
                for (j=0; j < Vars->Coord_Bin[i+1]; j++) /* search min of signal */
                  if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j] > 0)) XY = Vars->Mon2D_p[i][j];
                if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
              } /* if */

              for (j=0; j < Vars->Coord_Bin[i+1]; j++)
              {
                p1m[j] = Vars->Mon2D_p[i][j];
                p2m[j] = Vars->Mon2D_p2[i][j];
                if (Vars->Flag_signal != DEFS->COORD_P && Vars->Mon2D_N[i][j] > 0)
                { /* normalize mean signal to the number of events */
                  p1m[j] /= Vars->Mon2D_N[i][j];
                  p2m[j] /= Vars->Mon2D_N[i][j]*Vars->Mon2D_N[i][j];
                }
                if (Vars->Flag_log != 0)
                {
                  if ((p1m[j] > 0) && (p2m[j] > 0))
                  {
                    p2m[j] /= p1m[j]*p1m[j];
                    p1m[j] = log(p1m[j])/log(10);
                  }
                  else
                  {
                    p1m[j] = XY;
                    p2m[j] = 0;
                  }
                }
              } /* for */
              detector = mcdetector_out_1D(
                label,
                Vars->Coord_Label[i+1],
                Vars->Coord_Label[0],
                Vars->Coord_Var[i+1],
                min1d, max1d,
                Vars->Coord_Bin[i+1],
                Vars->Mon2D_N[i],p1m,p2m,
                fname, Vars->compcurname, Vars->compcurpos);

            } /* else */
            /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
            if (p1m != NULL) free(p1m); p1m=NULL;
            if (p2m != NULL) free(p2m); p2m=NULL;
            */
          } else { /* 0d monitor */
            detector = mcdetector_out_0D(label, Vars->Mon2D_p[i][0], Vars->Mon2D_p2[i][0], Vars->Mon2D_N[i][0], Vars->compcurname, Vars->compcurpos);
          }


        } /* for */
      } /* if 1D */
      else
      if (Vars->Coord_NumberNoPixel == 2)  /* 2D: DETECTOR_OUT_2D */
      {
        strcpy(fname,Vars->Mon_File);

        p0m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p1m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        p2m = (double *)malloc(Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
        if (p2m == NULL)
        {
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s cannot allocate memory for 2D array (%li). Skipping.\n", Vars->compcurname, 3*Vars->Coord_Bin[1]*Vars->Coord_Bin[2]*sizeof(double));
          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          */
        }
        else
        {
          if (Vars->Flag_log != 0)
          {
            XY = FLT_MAX;
            for (i= 0; i < Vars->Coord_Bin[1]; i++)
              for (j= 0; j < Vars->Coord_Bin[2]; j++) /* search min of signal */
                if ((XY > Vars->Mon2D_p[i][j]) && (Vars->Mon2D_p[i][j]>0)) XY = Vars->Mon2D_p[i][j];
            if (XY <= 0) XY = -log(FLT_MAX)/log(10); else XY = log(XY)/log(10)-1;
          }
          for (i= 0; i < Vars->Coord_Bin[1]; i++)
          {
            for (j= 0; j < Vars->Coord_Bin[2]; j++)
            {
              long index;
              index = j + i*Vars->Coord_Bin[2];
              p0m[index] = Vars->Mon2D_N[i][j];
              p1m[index] = Vars->Mon2D_p[i][j];
              p2m[index] = Vars->Mon2D_p2[i][j];
              if (Vars->Flag_signal != DEFS->COORD_P && p0m[index] > 0)
              {
                  p1m[index] /= p0m[index];
                  p2m[index] /= p0m[index]*p0m[index];
              }

              if (Vars->Flag_log != 0)
              {
                if ((p1m[index] > 0) && (p2m[index] > 0))
                {
                  p2m[index] /= (p1m[index]*p1m[index]);
                  p1m[index] = log(p1m[index])/log(10);

                }
                else
                {
                  p1m[index] = XY;
                  p2m[index] = 0;
                }
              }
            }
          }
          if (strchr(Vars->Mon_File,'.') == NULL)
          { strcat(fname, "."); strcat(fname, Vars->Coord_Var[1]);
              strcat(fname, "_"); strcat(fname, Vars->Coord_Var[2]); }
          if (Vars->Flag_Verbose) printf("Monitor_nD: %s write monitor file %s 2D (%lix%li).\n", Vars->compcurname, fname, Vars->Coord_Bin[1], Vars->Coord_Bin[2]);

          min1d = Vars->Coord_Min[1];
          max1d = Vars->Coord_Max[1];
          if (min1d == max1d) max1d = min1d+1e-6;
          min2d = Vars->Coord_Min[2];
          max2d = Vars->Coord_Max[2];
          if (min2d == max2d) max2d = min2d+1e-6;
          strcpy(label, Vars->Monitor_Label);
          if (Vars->Coord_Bin[1]*Vars->Coord_Bin[2] > 1
           && Vars->Flag_signal == DEFS->COORD_P)
            strcat(label, " per bin");

          detector = mcdetector_out_2D(
            label,
            Vars->Coord_Label[1],
            Vars->Coord_Label[2],
            min1d, max1d,
            min2d, max2d,
            Vars->Coord_Bin[1],
            Vars->Coord_Bin[2],
            p0m,p1m,p2m,
            fname, Vars->compcurname, Vars->compcurpos);

          /* comment out 'free memory' lines to avoid loosing arrays if
               'detector' structure is used by other instrument parts
          if (p0m != NULL) free(p0m);
          if (p1m != NULL) free(p1m);
          if (p2m != NULL) free(p2m);
          */
        }
      }
      free(fname);
    }
    return(detector);
  } /* end Monitor_nD_Save */

/* ========================================================================= */
/* Monitor_nD_Finally: this routine is used to free memory                   */
/* ========================================================================= */

void Monitor_nD_Finally(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    int i;

    /* Now Free memory Mon2D.. */
    if ((Vars->Flag_Auto_Limits || Vars->Flag_List) && Vars->Coord_Number)
    { /* Dim : (Vars->Coord_Number+1)*Vars->Buffer_Block matrix (for p, dp) */
      if (Vars->Mon2D_Buffer != NULL) free(Vars->Mon2D_Buffer);
    }

    /* 1D and n1D case : Vars->Flag_Multiple */
    if (Vars->Flag_Multiple && Vars->Coord_Number)
    { /* Dim : Vars->Coord_Number*Vars->Coord_Bin[i] vectors */
      for (i= 0; i < Vars->Coord_Number; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }


    /* 2D case : Vars->Coord_Number==2 and !Vars->Flag_Multiple and !Vars->Flag_List */
    if ((Vars->Coord_NumberNoPixel == 2) && !Vars->Flag_Multiple)
    { /* Dim : Vars->Coord_Bin[1]*Vars->Coord_Bin[2] matrix */
      for (i= 0; i < Vars->Coord_Bin[1]; i++)
      {
        free(Vars->Mon2D_N[i]);
        free(Vars->Mon2D_p[i]);
        free(Vars->Mon2D_p2[i]);
      }
      free(Vars->Mon2D_N);
      free(Vars->Mon2D_p);
      free(Vars->Mon2D_p2);
    }
  } /* end Monitor_nD_Finally */

/* ========================================================================= */
/* Monitor_nD_McDisplay: this routine is used to display component           */
/* ========================================================================= */

void Monitor_nD_McDisplay(MonitornD_Defines_type *DEFS,
  MonitornD_Variables_type *Vars)
  {
    double radius, h;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    int    i;
    double hdiv_min=-180, hdiv_max=180, vdiv_min=-90, vdiv_max=90;
    char   restricted = 0;

    radius = Vars->Sphere_Radius;
    h = Vars->Cylinder_Height;
    xmin = Vars->mxmin;
    xmax = Vars->mxmax;
    ymin = Vars->mymin;
    ymax = Vars->mymax;
    zmin = Vars->mzmin;
    zmax = Vars->mzmax;

    /* determine if there are angular limits set at start (no auto) in coord_types
     * cylinder/banana: look for hdiv
     * sphere: look for angle, radius (->atan2(val,radius)), hdiv, vdiv
     * this activates a 'restricted' flag, to draw a region as blades on cylinder/sphere
     */
    for (i= 0; i <= Vars->Coord_Number; i++)
    {
      int Set_Vars_Coord_Type;
      Set_Vars_Coord_Type = (Vars->Coord_Type[i] & (DEFS->COORD_LOG-1));
      if (Set_Vars_Coord_Type == DEFS->COORD_HDIV || Set_Vars_Coord_Type == DEFS->COORD_THETA)
      { hdiv_min = Vars->Coord_Min[i]; hdiv_max = Vars->Coord_Max[i]; restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_VDIV || Set_Vars_Coord_Type == DEFS->COORD_PHI)
      { vdiv_min = Vars->Coord_Min[i]; vdiv_max = Vars->Coord_Max[i];restricted = 1;  }
      else if (Set_Vars_Coord_Type == DEFS->COORD_ANGLE)
      { hdiv_min = vdiv_min = Vars->Coord_Min[i];
        hdiv_max = vdiv_max = Vars->Coord_Max[i];
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_RADIUS)
      { double angle;
        angle = RAD2DEG*atan2(Vars->Coord_Max[i], radius);
        hdiv_min = vdiv_min = angle;
        hdiv_max = vdiv_max = angle;
        restricted = 1; }
      else if (Set_Vars_Coord_Type == DEFS->COORD_Y && abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)
      {
        vdiv_min = atan2(ymin,radius)*RAD2DEG;
        vdiv_max = atan2(ymax,radius)*RAD2DEG;
        restricted = 1;
      }
    }
    /* full sphere */
    if ((!restricted && (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE))
    || abs(Vars->Flag_Shape) == DEFS->SHAPE_PREVIOUS)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
      mcdis_circle("xz",0,0,0,radius);
      mcdis_circle("yz",0,0,0,radius);
    }
    /* banana/cylinder/sphere portion */
    else
    if (restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)
                    || (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE)))
    {
      int NH=24, NV=24;
      int ih, iv;
      double width, height;
      int issphere;
      issphere = (abs(Vars->Flag_Shape) == DEFS->SHAPE_SPHERE);
      width = (hdiv_max-hdiv_min)/NH;
      if (!issphere) NV=1; /* cylinder has vertical axis */
      else height= (vdiv_max-vdiv_min)/NV;
      
      /* check width and height of elements (sphere) to make sure the nb
         of plates remains limited */
      if (width < 10  && NH > 1) { width = 10;  NH=(hdiv_max-hdiv_min)/width; width=(hdiv_max-hdiv_min)/NH; }
      if (height < 10 && NV > 1) { height = 10; NV=(vdiv_max-vdiv_min)/height; height= (vdiv_max-vdiv_min)/NV; }
      
      mcdis_magnify("xyz");
      for(ih = 0; ih < NH; ih++)
        for(iv = 0; iv < NV; iv++)
        {
          double theta0, phi0, theta1, phi1;          /* angles in spherical coordinates */
          double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3; /* vertices at plate edges */
          phi0 = (hdiv_min+ width*ih-90)*DEG2RAD;        /* in xz plane */
          phi1 = (hdiv_min+ width*(ih+1)-90)*DEG2RAD;
          if (issphere)
          {
            theta0= (vdiv_min+height* iv + 90)   *DEG2RAD; /* in vertical plane */
            theta1= (vdiv_min+height*(iv+1) + 90)*DEG2RAD;
            
            y0 = -radius*cos(theta0);            /* z with Z vertical */
            y1 = -radius*cos(theta1);
            if (y0 < ymin) y0=ymin;
            if (y0 > ymax) y0=ymax;
            if (y1 < ymin) y1=ymin;
            if (y1 > ymax) y1=ymax;
          } else {
            y0 = ymin;
            y1 = ymax;
            theta0=theta1=90*DEG2RAD;
          }

          x0 = radius*sin(theta0)*cos(phi0); /* x with Z vertical */
          z0 =-radius*sin(theta0)*sin(phi0); /* y with Z vertical */
          x1 = radius*sin(theta1)*cos(phi0); 
          z1 =-radius*sin(theta1)*sin(phi0);
          x2 = radius*sin(theta1)*cos(phi1); 
          z2 =-radius*sin(theta1)*sin(phi1);
          x3 = radius*sin(theta0)*cos(phi1); 
          z3 =-radius*sin(theta0)*sin(phi1);
          y2 = y1; y3 = y0;

          mcdis_multiline(5,
            x0,y0,z0,
            x1,y1,z1,
            x2,y2,z2,
            x3,y3,z3,
            x0,y0,z0);
        }
      if (Vars->Flag_mantid) {
	/* First define the base pixel type */
	double dt, dy;
	dt = (Vars->Coord_Max[1]-Vars->Coord_Min[1])/Vars->Coord_Bin[1];
	dy = (Vars->Coord_Max[2]-Vars->Coord_Min[2])/Vars->Coord_Bin[2];
	printf("MANTID_BANANA_DET:  %g, %g, %g, %g, %g, %li, %li, %g\n", radius, 
	       Vars->Coord_Min[1],Vars->Coord_Max[1], Vars->Coord_Min[2],Vars->Coord_Max[2], Vars->Coord_Bin[1], Vars->Coord_Bin[2], Vars->Coord_Min[4]); 
      }
    }
    /* disk (circle) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_DISK)
    {
      mcdis_magnify("");
      mcdis_circle("xy",0,0,0,radius);
    }
    /* rectangle (square) */
    else
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_SQUARE)
    {
      mcdis_magnify("xy");
      mcdis_multiline(5, (double)xmin, (double)ymin, 0.0,
             (double)xmax, (double)ymin, 0.0,
             (double)xmax, (double)ymax, 0.0,
             (double)xmin, (double)ymax, 0.0,
             (double)xmin, (double)ymin, 0.0);
      
      if (Vars->Flag_mantid) {
	/* First define the base pixel type */
	double dx, dy;
	dx = (Vars->Coord_Max[1]-Vars->Coord_Min[1])/Vars->Coord_Bin[1];
	dy = (Vars->Coord_Max[2]-Vars->Coord_Min[2])/Vars->Coord_Bin[2];
	printf("MANTID_RECTANGULAR_DET:  %g, %g, %g, %g, %li, %li, %g\n", 
	       Vars->Coord_Min[1],Vars->Coord_Max[1], Vars->Coord_Min[2],Vars->Coord_Max[2], Vars->Coord_Bin[1], Vars->Coord_Bin[2], Vars->Coord_Min[4]);
      }
    }
    /* full cylinder/banana */
    else
    if (!restricted && ((abs(Vars->Flag_Shape) == DEFS->SHAPE_CYLIND) || (abs(Vars->Flag_Shape) == DEFS->SHAPE_BANANA)))
    {
      mcdis_magnify("xyz");
      mcdis_circle("xz", 0,  h/2.0, 0, radius);
      mcdis_circle("xz", 0, -h/2.0, 0, radius);
      mcdis_line(-radius, -h/2.0, 0, -radius, +h/2.0, 0);
      mcdis_line(+radius, -h/2.0, 0, +radius, +h/2.0, 0);
      mcdis_line(0, -h/2.0, -radius, 0, +h/2.0, -radius);
      mcdis_line(0, -h/2.0, +radius, 0, +h/2.0, +radius);
    }
    else
    /* box */
    if (abs(Vars->Flag_Shape) == DEFS->SHAPE_BOX)
    {
      mcdis_magnify("xyz");
      mcdis_multiline(5, xmin, ymin, zmin,
                   xmax, ymin, zmin,
                   xmax, ymax, zmin,
                   xmin, ymax, zmin,
                   xmin, ymin, zmin);
      mcdis_multiline(5, xmin, ymin, zmax,
                   xmax, ymin, zmax,
                   xmax, ymax, zmax,
                   xmin, ymax, zmax,
                   xmin, ymin, zmax);
      mcdis_line(xmin, ymin, zmin, xmin, ymin, zmax);
      mcdis_line(xmax, ymin, zmin, xmax, ymin, zmax);
      mcdis_line(xmin, ymax, zmin, xmin, ymax, zmax);
      mcdis_line(xmax, ymax, zmin, xmax, ymax, zmax);
    }
  } /* end Monitor_nD_McDisplay */

/* end of monitor_nd-lib.c */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
    char    quiet;   /*output level for messages to the console 0: print all messages, 1:only print some/including errors, 2: never print anything.*/
  } t_Table;

/*maximum number of rows to rebin a table = 1M*/
enum { mcread_table_rebin_maxsize = 1000000 };

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
#pragma acc routine
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
#pragma acc routine
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
#pragma acc routine
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
#pragma acc routine
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
#pragma acc routine
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
#pragma acc routine
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);


#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value should just be used as
     * if the table had been read from disk. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we should proceed with freeing the memory
     * associated with the table item - otherwise only decrement the reference counter since there are more references
     * that may need it.*/

    static t_Read_table_file_item read_table_file_list[1024];  
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref = ((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item - 
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data && 
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found and no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found and the reference counter is 1.
                         * This means we should garbage collect. Move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            /* item not found, and so should be garbage collected. This could be the case if freeing a
             * Table that has been constructed from code - not read from file. Return 0x1 to flag it for
             * collection.*/
            return (void *) 0x1 ;
    }
    /* If we arrive here, nothing worked, return NULL */
    return NULL;
}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None. 
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    return Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;
    
    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);
    
    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && instrument_source[0] != '\0' && strlen(instrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(instrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - instrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, instrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && instrument_exe[0] != '\0' && strlen(instrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(instrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - instrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, instrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        // fprintf(stderr, "Warning: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);
    
    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else                   
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        *Table=*tab_p;
        MPI_MASTER(
            if(Table->quiet<1)
              printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
            );
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
          if(Table->quiet<1)
            printf("Opening input file '%s' (Table_Read_Offset)\n", path);
          );
    }
    
    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    
    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    
    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    
    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
          if(Table->quiet<1)
            printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }
    
    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    
    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      if(!(Table->quiet>1))
        fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      if(!(Table->quiet>1))
        fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  count_invalid       = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }
        if (strstr(line, "***"))
        {
          count_invalid++;
          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array*1.5;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename[0] != '\0' ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    if (count_invalid)
    {
      fprintf(stderr,"Warning: Read_Table :%s %s Data has %i invalid lines (*****). Ignored.\n",
      (Table->filename[0] != '\0' ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_invalid);
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data 
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      /*return early if the rebinned table will become too large*/
      if (Length_Table > mcread_table_rebin_maxsize){
        fprintf(stderr,"WARNING: (Table_Rebin): Rebinning table from %s would exceed 1M rows. Skipping.\n", Table->filename); 
        return(Table->rows*Table->columns);
      }
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step; 
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index-1,0);
    X2 = Table_Index(Table,Index  ,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index <= Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1, j);
  Y2 = Table_Index(Table,Index  , j);

#ifdef OPENACC
#define strcmp(a,b) str_comp(a,b)
#endif

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

#ifdef OPENACC
#ifdef strcmp
#undef strcmp
#endif
#endif

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

#ifdef OPENACC
#define strcmp(a,b) str_comp(a,b)
#endif

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
#ifdef OPENACC
#ifdef strcmp
#undef strcmp
#endif
#endif
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table. First Call Table_File_list_gc. If this returns
*   non-zero it means there are more refernces to the table, and so the table
*   should not bee freed.
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    } 
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
        Table.filename[0] != '\0' ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      if(Table->quiet<2)
        fprintf(stderr,"Error: allocating %ld double elements."
            "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* first allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      
      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /*if the block is empty - don't store it*/
      if (nelements>0){
          /* if t_Table array is not long enough, expand and realocate */
          if (block_number >= allocated-1) {
              allocated += 256;
              Table_Array = (t_Table *)realloc(Table_Array,
                      allocated*sizeof(t_Table));
              if (!Table_Array) {
                  fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
                          allocated*sizeof(t_Table));
                  *blocks = 0;
                  return (NULL);
              }
          }
          /* store it into t_Table array */
          //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
          Table_Array[block_number-1] = Table;
      }
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index;
    if (!Table) return;
    for (index=0;index < Table[0].array_length; index++){
            Table_Free(&Table[index]);
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interoff.h
*
* %Identification
* Written by: Reynald Arnerin
* Date:    Jun 12, 2008
* Release:
* Version:
*
* Object File Format intersection header for McStas. Requires the qsort function.
*
* Such files may be obtained with e.g.
*   qhull < points.xyz Qx Qv Tv o > points.off
* where points.xyz has format:
*   3
*   <nb_points>
*   <x> <y> <z>
*   ...
* The resulting file should have its first line being changed from '3' into 'OFF'.
* It can then be displayed with geomview.
* A similar, but somewhat older solution is to use 'powercrust' with e.g.
*   powercrust -i points.xyz
* which will generate a 'pc.off' file to be renamed as suited.
*
*******************************************************************************/

#ifndef INTEROFF_LIB_H
#define INTEROFF_LIB_H "$Revision$"

#ifndef OFF_EPSILON
#define OFF_EPSILON 1e-13
#endif

#ifndef OFF_INTERSECT_MAX
#ifdef OPENACC
#define OFF_INTERSECT_MAX 100
#else
#define OFF_INTERSECT_MAX 1024
#endif
#endif

//#include <float.h>

#define N_VERTEX_DISPLAYED    200000

typedef struct intersection {
	MCNUM time;  	  //time of the intersection
	Coords v;	      //intersection point
	Coords normal;  //normal vector of the surface intersected
	short in_out;	  //1 if the ray enters the volume, -1 otherwise
	short edge;	    //1 if the intersection is on the boundary of the polygon, and error is possible
	unsigned long index; // index of the face
} intersection;

typedef struct polygon {
  MCNUM* p;       //vertices of the polygon in adjacent order, this way : x1 | y1 | z1 | x2 | y2 | z2 ...
  int npol;       //number of vertices
  #pragma acc shape(p[0:npol]) init_needed(npol)
  Coords normal;
  double D;
} polygon;

typedef struct off_struct {
    long vtxSize;
    long polySize;
    long faceSize;
    Coords* vtxArray;
    #pragma acc shape(vtxArray[0:vtxSize]) init_needed(vtxSize)
    Coords* normalArray;
    #pragma acc shape(vtxArray[0:faceSize]) init_needed(faceSize)
    unsigned long* faceArray;
    #pragma acc shape(vtxArray[0:faceSize][0:polySize]) init_needed(faceSize,polySize)
    double* DArray;
    #pragma acc shape(vtxArray[0:polySize]) init_needed(polySize)
    char *filename;
    int mantidflag;
    long mantidoffset;
    intersection intersects[OFF_INTERSECT_MAX]; // After a call to off_intersect_all contains the list of intersections.
    int nextintersect;                 // 'Next' intersection (first t>0) solution after call to off_intersect_all
    int numintersect;               // Number of intersections after call to off_intersect_all
} off_struct;

/*******************************************************************************
* long off_init(  char *offfile, double xwidth, double yheight, double zdepth, off_struct* data)
* ACTION: read an OFF file, optionally center object and rescale, initialize OFF data structure
* INPUT: 'offfile' OFF file to read
*        'xwidth,yheight,zdepth' if given as non-zero, apply bounding box.
*           Specifying only one of these will also use the same ratio on all axes
*        'notcenter' center the object to the (0,0,0) position in local frame when set to zero
* RETURN: number of polyhedra and 'data' OFF structure
*******************************************************************************/
long off_init(  char *offfile, double xwidth, double yheight, double zdepth,
                int notcenter, off_struct* data);

/*******************************************************************************
* int off_intersect_all(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct *data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         ax, ay, az are the local acceleration vector
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*         data is the full OFF structure, including a list intersection type
*******************************************************************************/
#pragma acc routine
int off_intersect_all(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct *data );

/*******************************************************************************
* int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         ax, ay, az are the local acceleration vector
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
#pragma acc routine
int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct data );

/*****************************************************************************
* int off_intersectx(double* l0, double* l3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double kx, double ky, double kz,
     off_struct data )
* ACTION: computes intersection of an xray trajectory with an object.
* INPUT:  x,y,z and kx,ky,kz, are spatial coordinates and wavevector of the x-ray
*         respectively. data points to the OFF data structure.
* RETURN: the number of polyhedra the trajectory intersects
*         l0 and l3 are the smallest incoming and outgoing intersection lengths
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
#pragma acc routine
int off_x_intersect(double *l0,double *l3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double kx, double ky, double kz,
     off_struct data );

/*******************************************************************************
* void off_display(off_struct data)
* ACTION: display up to N_VERTEX_DISPLAYED points from the object
*******************************************************************************/
void off_display(off_struct);

/*******************************************************************************
void p_to_quadratic(double eq[], Coords acc,
                    Coords pos, Coords vel,
                    double* teq)
* ACTION: define the quadratic for the intersection of a parabola with a plane
* INPUT: 'eq' plane equation
*        'acc' acceleration vector
*        'vel' velocity of the particle
*        'pos' position of the particle
*         equation of plane A * x + B * y + C * z - D = 0
*         eq[0] = (C*az)/2+(B*ay)/2+(A*ax)/2
*         eq[1] = C*vz+B*vy+A*vx
*         eq[2] = C*z0+B*y0+A*x0-D
* RETURN: equation of parabola: teq(0) * t^2 + teq(1) * t + teq(2)
*******************************************************************************/
void p_to_quadratic(Coords norm, MCNUM d, Coords acc, Coords pos, Coords vel,
		    double* teq);

/*******************************************************************************
int quadraticSolve(double eq[], double* x1, double* x2);
* ACTION: solves the quadratic for the roots x1 and x2 
*         eq[0] * t^2 + eq[1] * t + eq[2] = 0
* INPUT: 'eq' the coefficients of the parabola
* RETURN: roots x1 and x2 and the number of solutions
*******************************************************************************/
int quadraticSolve(double* eq, double* x1, double* x2);

#endif

/* end of interoff-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2008, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/interoff-lib.c
*
* %Identification
* Written by: Reynald Arnerin
* Date:    Jun 12, 2008
* Origin: ILL
* Release: $Revision$
* Version: McStas X.Y
*
* Object File Format intersection library for McStas. Requires the qsort function.
*
* Such files may be obtained with e.g.
*   qhull < points.xyz Qx Qv Tv o > points.off
* where points.xyz has format (it supports comments):
*   3
*   <nb_points>
*   <x> <y> <z>
*   ...
* The resulting file should have its first line being changed from '3' into 'OFF'.
* It can then be displayed with geomview.
* A similar, but somewhat older solution is to use 'powercrust' with e.g.
*   powercrust -i points.xyz
* which will generate a 'pc.off' file to be renamed as suited.
*
*******************************************************************************/

#ifndef INTEROFF_LIB_H
#include "interoff-lib.h"
#endif

#pragma acc routine
double off_F(double x, double y,double z,double A,double B,double C,double D) {
  return ( A*x + B*y + C*z + D );
}

#pragma acc routine
char off_sign(double a) {
  if (a<0)       return(-1);
  else if (a==0) return(0);
  else           return(1);
}

// off_normal ******************************************************************
//gives the normal vector of p
#pragma acc routine
void off_normal(Coords* n, polygon p)
{
  //using Newell method
  int i=0,j=0;
  n->x=0;n->y=0;n->z=0;
  for (i = 0, j = p.npol-1; i < p.npol; j = i++)
  {
    MCNUM x1=p.p[3*i],
          y1=p.p[3*i+1],
          z1=p.p[3*i+2];
    MCNUM x2=p.p[3*j],
          y2=p.p[3*j+1],
          z2=p.p[3*j+2];
    // n is the cross product of v1*v2
    n->x += (y1 - y2) * (z1 + z2);
    n->y += (z1 - z2) * (x1 + x2);
    n->z += (x1 - x2) * (y1 + y2);
  }
} /* off_normal */

// off_pnpoly ******************************************************************
//based on http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//return 0 if the vertex is out
//    1 if it is in
//   -1 if on the boundary
#pragma acc routine
int off_pnpoly(polygon p, Coords v)
{
  int i=0, c = 0;
  MCNUM minx=FLT_MAX,maxx=-FLT_MAX,miny=FLT_MAX,maxy=-FLT_MAX,minz=FLT_MAX,maxz=-FLT_MAX;
  MCNUM areax=0,areay=0,areaz=0;

  int pol2dx=0,pol2dy=1;          //2d restriction of the poly
  MCNUM x=v.x,y=v.y;

  /*areax: projected area with x-scratched = |v1_yz x v2_yz|, where v1=(x1-x0,0,z1-z0) & v2=(x2-x0,0,z2-z0).*/
  /* In principle, if polygon is triangle area should be scaled by 1/2, but this is irrelevant for finding the maximum area.*/
  /* Similarly for y and z scratched.*/
  areax=coords_len(coords_xp(
        coords_set(0,p.p[3*1+1]-p.p[0+1],p.p[3*1+2]-p.p[0+2]),
        coords_set(0,p.p[3*2+1]-p.p[0+1],p.p[3*2+2]-p.p[0+2])));
  areay=coords_len(coords_xp(
        coords_set(p.p[3*1+0]-p.p[0+0],0,p.p[3*1+2]-p.p[0+2]),
        coords_set(p.p[3*2+0]-p.p[0+0],0,p.p[3*2+2]-p.p[0+2])));
  areaz=coords_len(coords_xp(
        coords_set(p.p[3*1+0]-p.p[0+0],p.p[3*1+1]-p.p[0+1],0),
        coords_set(p.p[3*2+0]-p.p[0+0],p.p[3*2+1]-p.p[0+1],0)));

  if(areaz<areax){
    if(areax<areay){
      /*pick areay - i.e. scratch y*/
      pol2dy=2;
      y=v.z;
    }else{
      /*scratch x*/
      pol2dx=2;
      x=v.z;
    }
  }else if (areaz<areay){
    pol2dy=2;
    y=v.z;
  }

  //trace rays and test number of intersection
  int j;
  for (i = 0, j = p.npol-1; i < p.npol; j = i++) {
    if (((((p.p[3*i+pol2dy])<=y) && (y<(p.p[3*j+pol2dy]))) ||
         (((p.p[3*j+pol2dy])<=y) && (y<(p.p[3*i+pol2dy])))) &&
        (x < ( (p.p[3*j+pol2dx] - p.p[3*i+pol2dx]) * (y - p.p[3*i+pol2dy])
             / (p.p[3*j+pol2dy] - p.p[3*i+pol2dy]) + p.p[3*i+pol2dx]) ))
      c = !c;

    if (((fabs(p.p[3*i+pol2dy]-y)<=OFF_EPSILON) || ((fabs(p.p[3*j+pol2dy]-y)<=OFF_EPSILON))) &&
        fabs(x -((p.p[3*j+pol2dx] - p.p[3*i+pol2dx]) * (y - p.p[3*i+pol2dy])
          / (p.p[3*j+pol2dy] - p.p[3*i+pol2dy]) + p.p[3*i+pol2dx])) < OFF_EPSILON)
    {
      //the point lies on the edge
      c=-1;
      break;
    }
  }

  return c;
} /* off_pnpoly */

// off_intersectPoly ***********************************************************
//gives the intersection vertex between ray [a,b) and polygon p and its parametric value on (a b)
//based on http://geometryalgorithms.com/Archive/algorithm_0105/algorithm_0105.htm
#pragma acc routine
int off_intersectPoly(intersection *inter, Coords a, Coords b, polygon p)
{
  //direction vector of [a,b]
  Coords dir = {b.x-a.x, b.y-a.y, b.z-a.z};

  //the normal vector to the polygon
  Coords normale=p.normal;
  //off_normal(&normale, p); done at the init stage

  //direction vector from a to a vertex of the polygon
  Coords w0 = {a.x-p.p[0], a.y-p.p[1], a.z-p.p[2]};

  //scalar product
  MCNUM nw0  =-scalar_prod(normale.x,normale.y,normale.z,w0.x,w0.y,w0.z);
  MCNUM ndir = scalar_prod(normale.x,normale.y,normale.z,dir.x,dir.y,dir.z);
  inter->time = inter->edge = inter->in_out=0;
  inter->v = inter->normal = coords_set(0,0,1);

  if (fabs(ndir) < OFF_EPSILON)    // ray is parallel to polygon plane
  {
    if (nw0 == 0)              // ray lies in polygon plane (infinite number of solution)
      return 0;
    else return 0;             // ray disjoint from plane (no solution)
  }

  // get intersect point of ray with polygon plane
  inter->time = nw0 / ndir;            //parametric value the point on line (a,b)

  inter->v = coords_set(a.x + inter->time * dir.x,// intersect point of ray and plane
    a.y + inter->time * dir.y,
    a.z + inter->time * dir.z);

  int res=off_pnpoly(p,inter->v);

  inter->edge=(res==-1);
  if (ndir<0)
    inter->in_out=1;  //the negative dot product means we enter the surface
  else
    inter->in_out=-1;

  inter->normal=p.normal;

  return res;         //true if the intersection point lies inside the poly
} /* off_intersectPoly */


// off_getBlocksIndex **********************************************************
/*reads the indexes at the beginning of the off file as this :
line 1  OFF
line 2  nbVertex nbFaces nbEdges
*/
FILE *off_getBlocksIndex(char* filename, long* vtxSize, long* polySize )
{
  FILE* f = Open_File(filename,"r", NULL); /* from read_table-lib: FILE *Open_File(char *name, char *Mode, char *path) */
  if (!f) return (f);

  char line[CHAR_BUF_LENGTH];
  char *ret=0;
  *vtxSize = *polySize = 0;

  /* **************** start to read the file header */
  /* OFF file:
     'OFF' or '3'
   */

  ret=fgets(line,CHAR_BUF_LENGTH , f);// line 1 = "OFF"
  if (ret == NULL)
  {
    fprintf(stderr, "Error: Can not read 1st line in file %s (interoff/off_getBlocksIndex)\n", filename);
    exit(1);
  }
  if (strlen(line)>5)
  {
      fprintf(stderr,"Error: First line in %s is too long (=%lu). Possibly the line is not terminated by '\\n'.\n"
              "       The first line is required to be exactly 'OFF', '3' or 'ply'.\n",
              filename,(long unsigned)strlen(line));
      fclose(f);
      return(NULL);
  }

  if (strncmp(line,"OFF",3) && strncmp(line,"3",1) && strncmp(line,"ply",1))
  {
    fprintf(stderr, "Error: %s is probably not an OFF, NOFF or PLY file (interoff/off_getBlocksIndex).\n"
                    "       Requires first line to be 'OFF', '3' or 'ply'.\n",filename);
    fclose(f);
    return(NULL);
  }

  if (!strncmp(line,"OFF",3) || !strncmp(line,"3",1)) {
    do  /* OFF file: skip # comments which may be there */
    {
      ret=fgets(line,CHAR_BUF_LENGTH , f);
      if (ret == NULL)
      {
        fprintf(stderr, "Error: Can not read line in file %s (interoff/off_getBlocksIndex)\n", filename);
        exit(1);
      }
    } while (line[0]=='#');
    //line = nblines of vertex,faces and edges arrays
    sscanf(line,"%lu %lu",vtxSize,polySize);
  } else {
    do  /* PLY file: read all lines until find 'end_header'
           and locate 'element faces' and 'element vertex' */
    {
      ret=fgets(line,CHAR_BUF_LENGTH , f);
      if (ret == NULL)
      {
        fprintf(stderr, "Error: Can not read line in file %s (interoff/off_getBlocksIndex)\n", filename);
        exit(1);
      }
      if (!strncmp(line,"element face",12))
        sscanf(line,"element face %lu",polySize);
      else if (!strncmp(line,"element vertex",14))
        sscanf(line,"element vertex %lu",vtxSize);
      else if (!strncmp(line,"format binary",13))
        exit(fprintf(stderr,
          "Error: Can not read binary PLY file %s, only 'format ascii' (interoff/off_getBlocksIndex)\n%s\n",
          filename, line));
    } while (strncmp(line,"end_header",10));
  }

  /* The FILE is left opened ready to read 'vtxSize' vertices (vtxSize *3 numbers)
     and then polySize polygons (rows) */

  return(f);
} /* off_getBlocksIndex */

// off_init_planes *************************************************************
//gives the equations of 2 perpandicular planes of [ab]
#pragma acc routine
void off_init_planes(Coords a, Coords b,
  MCNUM* A1, MCNUM* C1, MCNUM* D1, MCNUM *A2, MCNUM* B2, MCNUM* C2, MCNUM* D2)
{
  //direction vector of [a b]
  Coords dir={b.x-a.x, b.y-a.y, b.z-a.z};

  //the plane parallel to the 'y' is computed with the normal vector of the projection of [ab] on plane 'xz'
  *A1= dir.z;
  *C1=-dir.x;
  if(*A1!=0 || *C1!=0)
    *D1=-(a.x)*(*A1)-(a.z)*(*C1);
  else
  {
    //the plane does not support the vector, take the one parallel to 'z''
    *A1=1;
    //B1=dir.x=0
    *D1=-(a.x);
  }
  //the plane parallel to the 'x' is computed with the normal vector of the projection of [ab] on plane 'yz'
  *B2= dir.z;
  *C2=-dir.y;
  *A2= 0;
  if (*B2==0 && *C2==0)
  {
    //the plane does not support the vector, take the one parallel to 'z'
    *B2=1;
    //B1=dir.x=0
    *D2=-(a.y);
  }
  else {
    if (dir.z==0)
    {
      //the planes are the same, take the one parallel to 'z'
      *A2= dir.y;
      *B2=-dir.x;
      *D2=-(a.x)*(*A2)-(a.y)*(*B2);
    }
    else
      *D2=-(a.y)**B2-(a.z)**C2;
  }
} /* off_init_planes */

// off_clip_3D_mod *************************************************************
#pragma acc routine
int off_clip_3D_mod(intersection* t, Coords a, Coords b,
  Coords* vtxArray, unsigned long vtxSize, unsigned long* faceArray,
  unsigned long faceSize, Coords* normalArray)
{
  MCNUM A1=0, C1=0, D1=0, A2=0, B2=0, C2=0, D2=0;      //perpendicular plane equations to [a,b]
  off_init_planes(a, b, &A1, &C1, &D1, &A2, &B2, &C2, &D2);

  int t_size=0;
  MCNUM popol[3*4]; /*3 dimensions and max 4 vertices to form a polygon*/
  unsigned long i=0,indPoly=0;

  //exploring the polygons :
  i=indPoly=0;
  while (i<faceSize)
  {
    polygon pol;
    pol.npol  = faceArray[i];                //nb vertex of polygon
    pol.p     = popol;
    pol.normal= coords_set(0,0,1);
    unsigned long indVertP1=faceArray[++i];  //polygon's first vertex index in vtxTable
    int j=1;
    /*check whether vertex is left or right of plane*/
    char sg0=off_sign(off_F(vtxArray[indVertP1].x,vtxArray[indVertP1].y,vtxArray[indVertP1].z,A1,0,C1,D1));
    while (j<pol.npol)
    {
      //polygon's j-th vertex index in vtxTable
      unsigned long indVertP2=faceArray[i+j];
      /*check whether vertex is left or right of plane*/
      char sg1=off_sign(off_F(vtxArray[indVertP2].x,vtxArray[indVertP2].y,vtxArray[indVertP2].z,A1,0,C1,D1));
      if (sg0!=sg1) //if the plane intersect the polygon
        break;

      ++j;
    }

    if (j<pol.npol)          //ok, let's test with the second plane
    {
      char sg1=off_sign(off_F(vtxArray[indVertP1].x,vtxArray[indVertP1].y,vtxArray[indVertP1].z,A2,B2,C2,D2));//tells if vertex is left or right of the plane

      j=1;
      while (j<pol.npol)
      {
        //unsigned long indVertPi=faceArray[i+j];  //polyg's j-th vertex index in vtxTable
        Coords vertPi=vtxArray[faceArray[i+j]];
        if (sg1!=off_sign(off_F(vertPi.x,vertPi.y,vertPi.z,A2,B2,C2,D2)))//if the plane intersect the polygon
          break;
        ++j;
      }
      if (j<pol.npol)
      {
#ifdef OFF_LEGACY
        if (t_size>OFF_INTERSECT_MAX)
        {
#ifndef OPENACC
          fprintf(stderr, "Warning: number of intersection exceeded (%d) (interoff-lib/off_clip_3D_mod)\n", OFF_INTERSECT_MAX);
#endif
            return (t_size);
        }
#endif
        //both planes intersect the polygon, let's find the intersection point
        //our polygon :
        int k;
        for (k=0; k<pol.npol; ++k)
        {
          Coords vertPk=vtxArray[faceArray[i+k]];
          pol.p[3*k]  =vertPk.x;
          pol.p[3*k+1]=vertPk.y;
          pol.p[3*k+2]=vertPk.z;
        }
        pol.normal=normalArray[indPoly];
        intersection x;
        if (off_intersectPoly(&x, a, b, pol))
        {
          x.index = indPoly;
#ifdef OFF_LEGACY
          t[t_size++]=x;
#else
	  /* Check against our 4 existing times, starting from [-FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX] */
	  /* Case 1, negative time? */
	  if (t_size < 4) t_size++;	  
	  if (x.time < 0) {
	    if (x.time > t[0].time) {
	      t[0]=x;
	    }
	  } else {
	    /* Case 2, positive time */
	    intersection xtmp;
	    if (x.time < t[3].time) {
	      t[3]=x;
	      if (t[3].time < t[2].time) {
		xtmp = t[2];
		t[2] = t[3];
		t[3] = xtmp;
	      }
	      if (t[2].time < t[1].time) {
		xtmp = t[1];
		t[1] = t[2];
		t[2] = xtmp;
	      }
	    } 
	  }
#endif
	}
      } /* if (j<pol.npol) */
    } /* if (j<pol.npol) */
    i += pol.npol;
    indPoly++;
  } /* while i<faceSize */
  return t_size;
} /* off_clip_3D_mod */

// off_clip_3D_mod_grav *************************************************************
/*******************************************************************************
version of off_clip_3D_mod_grav
*******************************************************************************/
#pragma acc routine seq
int off_clip_3D_mod_grav(intersection* t, Coords pos, Coords vel, Coords acc,
  Coords* vtxArray, unsigned long vtxSize, unsigned long* faceArray,
  unsigned long faceSize, Coords* normalArray, double* DArray)
{
  int t_size=0;
  MCNUM popol[3*CHAR_BUF_LENGTH];
  double plane_Eq [4];
  double quadratic [3];
  unsigned long i=0,indPoly=0;
  //exploring the polygons :
  i=indPoly=0;
  while (i<faceSize)
  {
    polygon pol;
    pol.npol  = faceArray[i];                //nb vertex of polygon
    pol.p     = popol;
    pol.normal= coords_set(0,0,1);
    unsigned long indVertP1=faceArray[++i];  //polygon's first vertex index in vtxTable
    
    if (t_size>CHAR_BUF_LENGTH)
      {
#ifndef OPENACC
	fprintf(stderr, "Warning: number of intersection exceeded (%d) (interoff-lib/off_clip_3D_mod)\n", CHAR_BUF_LENGTH);
#endif
	return (t_size);
      }
    //both planes intersect the polygon, let's find the intersection point
    //our polygon :
    int k;
    for (k=0; k<pol.npol; ++k)
      {
	Coords vertPk=vtxArray[faceArray[i+k]];
	pol.p[3*k]  =vertPk.x;
	pol.p[3*k+1]=vertPk.y;
	pol.p[3*k+2]=vertPk.z;
      }
    pol.normal=normalArray[indPoly];
    pol.D=DArray[indPoly];
    p_to_quadratic(pol.normal, pol.D, acc, pos, vel, quadratic);
    double x1, x2;
    int nsol = quadraticSolve(quadratic, &x1, &x2);

    if (nsol >= 1) {
      double time = 1.0e36;
      if (x1 < time && x1 > 0.0) {
	time = x1;
      }
      if (nsol == 2 && x2 < time && x2 > 0.0) {
	time = x2;
      }
      if (time != 1.0e36) {
	intersection inters;
	double t2 = time * time * 0.5;
	double tx = pos.x + time * vel.x;
	if (acc.x != 0.0) {
	  tx = tx + t2 * acc.x;
	}
	double ty = pos.y + time * vel.y;
	if (acc.y != 0.0) {
	  ty = ty + t2 * acc.y;
	}
	double tz = pos.z + time * vel.z;
	if (acc.z != 0.0) {
	  tz = tz + t2 * acc.z;
	}
	inters.v = coords_set(tx, ty, tz);
	Coords tvel = coords_set(vel.x + time * acc.x,
				 vel.y + time * acc.y,
				 vel.z + time * acc.z);
	inters.time = time;
	inters.normal = pol.normal;
	inters.index = indPoly;
	int res=off_pnpoly(pol,inters.v);
	if (res != 0) {
	  inters.edge=(res==-1);
	  MCNUM ndir = scalar_prod(pol.normal.x,pol.normal.y,pol.normal.z,tvel.x,tvel.y,tvel.z);
	  if (ndir<0) {
	    inters.in_out=1;  //the negative dot product means we enter the surface
	  } else {
	    inters.in_out=-1;
	  }
#ifdef OFF_LEGACY
          t[t_size++]=inters;
#else
    /* Check against our 4 existing times, starting from [-FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX] */
    /* Case 1, negative time? */
    if (t_size < 4) t_size++;
    if (inters.time < 0) {
      if (inters.time > t[0].time) {
        t[0]=inters;
      }
    } else {
      /* Case 2, positive time */
      intersection xtmp;
      if (inters.time < t[3].time) {
      t[3]=inters;
        if (t[3].time < t[2].time) {
    xtmp = t[2];
    t[2] = t[3];
    t[3] = xtmp;
        }
        if (t[2].time < t[1].time) {
    xtmp = t[1];
    t[1] = t[2];
    t[2] = xtmp;
        }
      }
    }
#endif
	}
      }
    }
    i += pol.npol;
    indPoly++;
  } /* while i<faceSize */
  return t_size;
} /* off_clip_3D_mod_grav */

// off_compare *****************************************************************
#pragma acc routine
int off_compare (void const *a, void const *b)
{
   intersection const *pa = a;
   intersection const *pb = b;

   return off_sign(pa->time - pb->time);
} /* off_compare */

// off_cleanDouble *************************************************************
//given an array of intersections throw those which appear several times
//returns 1 if there is a possibility of error
#pragma acc routine
int off_cleanDouble(intersection* t, int* t_size)
{
  int i=1;
  intersection prev=t[0];
  while (i<*t_size)
  {
    int j=i;
    //for each intersection with the same time
    while (j<*t_size && fabs(prev.time-t[j].time)<OFF_EPSILON)
    {
      //if the intersection is the exact same erase it
      if (prev.in_out==t[j].in_out)
      {
        int k;
        for (k=j+1; k<*t_size; ++k)
        {
          t[k-1]=t[k];
        }
        *t_size-=1;
      }
      else
        ++j;
    }
    prev=t[i];
    ++i;

  }
  return 1;
} /* off_cleanDouble */

// off_cleanInOut **************************************************************
//given an array of intesections throw those which enter and exit in the same time
//Meaning the ray passes very close to the volume
//returns 1 if there is a possibility of error
#pragma acc routine
int off_cleanInOut(intersection* t, int* t_size)
{
  int i=1;
  intersection prev=t[0];
  while (i<*t_size)
  {
    //if two intersection have the same time but one enters and the other exits erase both
    //(such intersections must be adjacent in the array : run off_cleanDouble before)
    if (fabs(prev.time-t[i].time)<OFF_EPSILON && prev.in_out!=t[i].in_out)
    {
      int j=0;
      for (j=i+1; j<*t_size; ++j)
      {
        t[j-2]=t[j];
      }
      *t_size-=2;
      prev=t[i-1];
    }
    else
    {
      prev=t[i];
      ++i;
    }
  }
  return (*t_size);
} /* off_cleanInOut */

/* PUBLIC functions ******************************************************** */

/*******************************************************************************
* long off_init(  char *offfile, double xwidth, double yheight, double zdepth, off_struct* data)
* ACTION: read an OFF file, optionally center object and rescale, initialize OFF data structure
* INPUT: 'offfile' OFF file to read
*        'xwidth,yheight,zdepth' if given as non-zero, apply bounding box.
*           Specifying only one of these will also use the same ratio on all axes
*        'notcenter' center the object to the (0,0,0) position in local frame when set to zero
* RETURN: number of polyhedra and 'data' OFF structure
*******************************************************************************/
long off_init(  char *offfile, double xwidth, double yheight, double zdepth,
                int notcenter, off_struct* data)
{
  // data to be initialized
  long    vtxSize =0, polySize=0, i=0, ret=0, faceSize=0;
  Coords* vtxArray        =NULL;
  Coords* normalArray     =NULL;
  double* DArray          =NULL;
  unsigned long* faceArray=NULL;
  FILE*   f               =NULL; /* the FILE with vertices and polygons */
  double minx=FLT_MAX,maxx=-FLT_MAX,miny=FLT_MAX,maxy=-FLT_MAX,minz=FLT_MAX,maxz=-FLT_MAX;

  // get the indexes
  if (!data) return(0);

  MPI_MASTER(
  printf("Loading geometry file (OFF/PLY): %s\n", offfile);
  );

  f=off_getBlocksIndex(offfile,&vtxSize,&polySize);
  if (!f) return(0);

  // read vertex table = [x y z | x y z | ...] =================================
  // now we read the vertices as 'vtxSize*3' numbers and store it in vtxArray
  MPI_MASTER(
  printf("  Number of vertices: %ld\n", vtxSize);
  );
  vtxArray   = malloc(vtxSize*sizeof(Coords));
  if (!vtxArray) return(0);
  i=0;
  while (i<vtxSize && ~feof(f))
  {
    double x,y,z;
    ret=fscanf(f, "%lg%lg%lg", &x,&y,&z);
    if (!ret) {
      // invalid line: we skip it (probably a comment)
      char line[CHAR_BUF_LENGTH];
      char *s=fgets(line, CHAR_BUF_LENGTH, f);
      continue;
    }
    if (ret != 3) {
      fprintf(stderr, "Error: can not read [xyz] coordinates for vertex %li in file %s (interoff/off_init). Read %li values.\n",
        i, offfile, ret);
      exit(2);
    }
    vtxArray[i].x=x;
    vtxArray[i].y=y;
    vtxArray[i].z=z;

    //bounding box
    if (vtxArray[i].x<minx) minx=vtxArray[i].x;
    if (vtxArray[i].x>maxx) maxx=vtxArray[i].x;
    if (vtxArray[i].y<miny) miny=vtxArray[i].y;
    if (vtxArray[i].y>maxy) maxy=vtxArray[i].y;
    if (vtxArray[i].z<minz) minz=vtxArray[i].z;
    if (vtxArray[i].z>maxz) maxz=vtxArray[i].z;
    i++; // inquire next vertex
  }

  // resizing and repositioning params
  double centerx=0, centery=0, centerz=0;
  if (!notcenter) {
    centerx=(minx+maxx)*0.5;
    centery=(miny+maxy)*0.5;
    centerz=(minz+maxz)*0.5;
  }

  double rangex=-minx+maxx,
         rangey=-miny+maxy,
         rangez=-minz+maxz;

  double ratiox=1,ratioy=1,ratioz=1;

  if (xwidth && rangex)
  {
    ratiox=xwidth/rangex;
    ratioy=ratiox;
    ratioz=ratiox;
  }

  if (yheight && rangey)
  {
    ratioy=yheight/rangey;
    if(!xwidth)  ratiox=ratioy;
    ratioz=ratioy;
  }

  if (zdepth && rangez)
  {
    ratioz=zdepth/rangez;
    if(!xwidth)  ratiox=ratioz;
    if(!yheight) ratioy=ratioz;
  }

  rangex *= ratiox;
  rangey *= ratioy;
  rangez *= ratioz;

  //center and resize the object
  for (i=0; i<vtxSize; ++i)
  {
    vtxArray[i].x=(vtxArray[i].x-centerx)*ratiox+(!notcenter ? 0 : centerx);
    vtxArray[i].y=(vtxArray[i].y-centery)*ratioy+(!notcenter ? 0 : centery);
    vtxArray[i].z=(vtxArray[i].z-centerz)*ratioz+(!notcenter ? 0 : centerz);
  }

  // read face table = [nbvertex v1 v2 vn | nbvertex v1 v2 vn ...] =============
  MPI_MASTER(
  printf("  Number of polygons: %ld\n", polySize);
  );
  normalArray= malloc(polySize*sizeof(Coords));
  faceArray  = malloc(polySize*10*sizeof(unsigned long)); // we assume polygons have less than 9 vertices
  DArray     = malloc(polySize*sizeof(double));
  if (!normalArray || !faceArray || !DArray) return(0);

  // fill faces
  faceSize=0;
  i=0;
  while (i<polySize && ~feof(f)) {
    int  nbVertex=0, j=0;
    // read the length of this polygon
    ret=fscanf(f, "%d", &nbVertex);
    if (!ret) {
      // invalid line: we skip it (probably a comment)
      char line[CHAR_BUF_LENGTH];
      char *s=fgets(line, CHAR_BUF_LENGTH, f);
      continue;
    }
    if (ret != 1) {
      fprintf(stderr, "Error: can not read polygon %li length in file %s (interoff/off_init)\n",
        i, offfile);
      exit(3);
    }
    if (faceSize > polySize*10) {
      fprintf(stderr, "Error: %li exceeded allocated polygon array[%li] in file %s (interoff/off_init)\n",
        faceSize, polySize*10, offfile);
    }
    faceArray[faceSize++] = nbVertex; // length of the polygon/face
    // then read the vertex ID's
    for (j=0; j<nbVertex; j++) {
      double vtx=0;
      ret=fscanf(f, "%lg", &vtx);
      faceArray[faceSize++] = vtx;   // add vertices index after length of polygon
    }
    i++;
  }

  // precomputes normals
  long indNormal=0;//index in polyArray
  i=0;    //index in faceArray
  while (i<faceSize)
  {
    int    nbVertex=faceArray[i];//nb of vertices of this polygon
    double vertices[3*nbVertex];
    int j;

    for (j=0; j<nbVertex; ++j)
    {
      unsigned long indVertPj=faceArray[i+j+1];
      vertices[3*j]  =vtxArray[indVertPj].x;
      vertices[3*j+1]=vtxArray[indVertPj].y;
      vertices[3*j+2]=vtxArray[indVertPj].z;
    }

    polygon p;
    p.p   =vertices;
    p.npol=nbVertex;
    off_normal(&(p.normal),p);

    normalArray[indNormal]=p.normal;
    p.D = scalar_prod(p.normal.x,p.normal.y,p.normal.z,
		      vertices[0],vertices[1],vertices[2]);
    DArray[indNormal]=p.D;

    i += nbVertex+1;
    indNormal++;

  }

  MPI_MASTER(
  if (ratiox!=ratioy || ratiox!=ratioz || ratioy!=ratioz)
    printf("Warning: Aspect ratio of the geometry %s was modified.\n"
           "         If you want to keep the original proportions, specifiy only one of the dimensions.\n",
           offfile);
  if ( xwidth==0 && yheight==0 && zdepth==0 ) {
    printf("Warning: Neither xwidth, yheight or zdepth are defined.\n"
	   "           The file-defined (non-scaled) geometry the OFF geometry %s will be applied!\n",
           offfile);
  }
  printf("  Bounding box dimensions for geometry %s:\n", offfile);
  printf("    Length=%f (%.3f%%)\n", rangex, ratiox*100);
  printf("    Width= %f (%.3f%%)\n", rangey, ratioy*100);
  printf("    Depth= %f (%.3f%%)\n", rangez, ratioz*100);
  );

  data->vtxArray   = vtxArray;
  data->normalArray= normalArray;
  data->DArray     = DArray;
  data->faceArray  = faceArray;
  data->vtxSize    = vtxSize;
  data->polySize   = polySize;
  data->faceSize   = faceSize;
  data->filename   = offfile;
  #ifdef OPENACC
  acc_attach((void *)&vtxArray);
  acc_attach((void *)&normalArray);
  acc_attach((void *)&faceArray);
  #endif

  return(polySize);
} /* off_init */

#pragma acc routine
int Min_int(int x, int y) {
  return (x<y)? x :y;
}

 
#pragma acc routine
void merge(intersection *arr, int l, int m, int r)
{
int i, j, k;
int n1 = m - l + 1;
int n2 =  r - m;

/* create temp arrays */
intersection *L, *R;
 L = (intersection *)malloc(sizeof(intersection) * n1);
 R = (intersection *)malloc(sizeof(intersection) * n2);
/* Copy data to temp arrays L[] and R[] */
 #pragma acc loop independent
for (i = 0; i < n1; i++)
    L[i] = arr[l + i];
 #pragma acc loop independent
for (j = 0; j < n2; j++)
    R[j] = arr[m + 1+ j];

/* Merge the temp arrays back into arr[l..r]*/
i = 0;
j = 0;
k = l;

while (i < n1 && j < n2)
{
    if (L[i].time <= R[j].time)
    {
        arr[k] = L[i];
        i++;
    }
    else
    {
        arr[k] = R[j];
        j++;
    }
    k++;
}

/* Copy the remaining elements of L[], if there are any */

while (i < n1)
{
    arr[k] = L[i];
    i++;
    k++;
}

/* Copy the remaining elements of R[], if there are any */
while (j < n2)
{
    arr[k] = R[j];
    j++;
    k++;
}
free(L);
free(R);
}


#ifdef USE_OFF
#pragma acc routine
void gpusort(intersection *arr, int size)
{
  int curr_size;  // For current size of subarrays to be merged
  // curr_size varies from 1 to n/2
  int left_start; // For picking starting index of left subarray
  // to be merged
  // pcopying (R[0:n2])
  {
    for (curr_size=1; curr_size<=size-1; curr_size = 2*curr_size)
      {
	// Pick starting point of different subarrays of current size
	for (left_start=0; left_start<size-1; left_start += 2*curr_size)
	  {
	    // Find ending point of left subarray. mid+1 is starting
	    // point of right
	    int mid = left_start + curr_size - 1;

	    int right_end = Min_int(left_start + 2*curr_size - 1, size-1);

	    // Merge Subarrays arr[left_start...mid] & arr[mid+1...right_end]
	    if (mid < right_end) merge(arr, left_start, mid, right_end);
	  }
      }
  }
}
#endif

/*******************************************************************************
void p_to_quadratic(double eq[], Coords acc,
                    Coords pos, Coords vel,
                    double* teq)
* ACTION: define the quadratic for the intersection of a parabola with a plane
* INPUT: 'eq' plane equation
*        'acc' acceleration vector
*        'vel' velocity of the particle
*        'pos' position of the particle
*         equation of plane A * x + B * y + C * z - D = 0
*         eq[0] = (C*az)/2+(B*ay)/2+(A*ax)/2
*         eq[1] = C*vz+B*vy+A*vx
*         eq[2] = C*z0+B*y0+A*x0-D
* RETURN: equation of parabola: teq(0) * t^2 + teq(1) * t + teq(2)
*******************************************************************************/
void p_to_quadratic(Coords norm, MCNUM d, Coords acc, Coords pos, Coords vel,
		    double* teq)
{
  teq[0] = scalar_prod(norm.x, norm.y, norm.z, acc.x, acc.y, acc.z) * 0.5;
  teq[1] = scalar_prod(norm.x, norm.y, norm.z, vel.x, vel.y, vel.z);
  teq[2] = scalar_prod(norm.x, norm.y, norm.z, pos.x, pos.y, pos.z) - d;
  return;
}

/*******************************************************************************
int quadraticSolve(double eq[], double* x1, double* x2);
* ACTION: solves the quadratic for the roots x1 and x2 
*         eq[0] * t^2 + eq[1] * t + eq[2] = 0
* INPUT: 'eq' the coefficients of the parabola
* RETURN: roots x1 and x2 and the number of solutions
*******************************************************************************/
int quadraticSolve(double* eq, double* x1, double* x2)
{
  if (eq[0] == 0.0) { // This is a linear equation
    if (eq[1] != 0.0) { // one solution
      *x1 = -eq[2]/eq[1];
      *x2 = 1.0e36;
      return 1;
    }else { // no solutions, 1.0e36 will be ignored.
      *x1 = 1.0e36;
      *x2 = 1.0e36;
      return 0;
    }
  }
  double delta = eq[1]*eq[1]-4.0*eq[0]*eq[2];
  if (delta < 0.0) { // no solutions, both are imaginary
    *x1 = 1.0e36;
    *x2 = 1.0e36;
    return 0;
  }
  double s = 1.0;
  if (eq[1] < 0) {
    s = -1.0;
  }
  *x1 = (-eq[1] - s * sqrt(delta))/(2.0*eq[0]);
  if (eq[0] != 0.0) { //two solutions
    *x2 = eq[2]/(eq[0]*(*x1));
    return 2;
  } else { //one solution
    *x2 = 1.0e36;
    return 1;
  }
}

/*******************************************************************************
* int off_intersect_all(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct *data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         ax, ay, az are the local acceleration vector
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*         data is the full OFF structure, including a list intersection type
*******************************************************************************/
int off_intersect_all(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct *data )
{

    int t_size = 0;
#ifdef OFF_LEGACY

    if(mcgravitation) {
      Coords pos={ x,  y,  z};
      Coords vel={vx, vy, vz};
      Coords acc={ax, ay, az};
      t_size=off_clip_3D_mod_grav(data->intersects, pos, vel, acc,
				  data->vtxArray, data->vtxSize, data->faceArray,
				  data->faceSize, data->normalArray, data->DArray );
    } else {
    ///////////////////////////////////
    // non-grav
      Coords A={x, y, z};
      Coords B={x+vx, y+vy, z+vz};
      t_size=off_clip_3D_mod(data->intersects, A, B,
			     data->vtxArray, data->vtxSize, data->faceArray,
			     data->faceSize, data->normalArray );
    }
    #ifndef OPENACC
    qsort(data->intersects, t_size, sizeof(intersection),  off_compare);
    #else
    #ifdef USE_OFF
    gpusort(data->intersects, t_size);
    #endif
    #endif
    off_cleanDouble(data->intersects, &t_size);
    off_cleanInOut(data->intersects,  &t_size);

    /*find intersections "closest" to 0 (favouring positive ones)*/
    if(t_size>0){
      int i=0;
      if(t_size>1) {
        for (i=1; i < t_size-1; i++){
          if (data->intersects[i-1].time > 0 && data->intersects[i].time > 0)
            break;
        }

	data->nextintersect=i-1;
	data->numintersect=t_size;

        if (t0) *t0 = data->intersects[i-1].time;
        if (n0) *n0 = data->intersects[i-1].normal;
        if (t3) *t3 = data->intersects[i].time;
        if (n3) *n3 = data->intersects[i].normal;
      } else {
        if (t0) *t0 = data->intersects[0].time;
	      if (n0) *n0 = data->intersects[0].normal;
      }
      /* should also return t[0].index and t[i].index as polygon ID */
      data->nextintersect=(data->intersects[data->nextintersect]).index;
      return t_size;
    }
#else
    intersection intersect4[4];
    intersect4[0].time=-FLT_MAX;
    intersect4[1].time=FLT_MAX;
    intersect4[2].time=FLT_MAX;
    intersect4[3].time=FLT_MAX;
    if(mcgravitation) {
      Coords pos={ x,  y,  z};
      Coords vel={vx, vy, vz};
      Coords acc={ax, ay, az};
      t_size=off_clip_3D_mod_grav(intersect4, pos, vel, acc,
				  data->vtxArray, data->vtxSize, data->faceArray,
				  data->faceSize, data->normalArray, data->DArray);
    } else {
    ///////////////////////////////////
    // non-grav
      Coords A={x, y, z};
      Coords B={x+vx, y+vy, z+vz};
      t_size=off_clip_3D_mod(intersect4, A, B,
	  data->vtxArray, data->vtxSize, data->faceArray, data->faceSize, data->normalArray );
    }
    if(t_size>0){
      int i=0;
      if (intersect4[0].time == -FLT_MAX) i=1;
      data->numintersect=t_size;
      if (t0) *t0 = intersect4[i].time;
      if (n0) *n0 = intersect4[i].normal;
      if (t3) *t3 = intersect4[i+1].time;
      if (n3) *n3 = intersect4[i+1].normal;

      if (intersect4[1].time == FLT_MAX)
      {
        if (t3) *t3 = 0.0;
      }

      /* should also return t[0].index and t[i].index as polygon ID */
      data->nextintersect=(int)intersect4[i].index;
      return t_size;
    }
#endif
    return 0;
} /* off_intersect */

/*******************************************************************************
* int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double vx, double vy, double vz,
     off_struct data )
* ACTION: computes intersection of neutron trajectory with an object.
* INPUT:  x,y,z and vx,vy,vz are the position and velocity of the neutron
*         data points to the OFF data structure
* RETURN: the number of polyhedra which trajectory intersects
*         t0 and t3 are the smallest incoming and outgoing intersection times
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_intersect(double* t0, double* t3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double vx, double vy, double vz,
     double ax, double ay, double az,
     off_struct data )
{
  return off_intersect_all(t0, t3, n0, n3, x, y, z, vx, vy, vz, ax, ay, az, &data );
} /* off_intersect */

/*****************************************************************************
* int off_x_intersect(double* l0, double* l3,
     Coords *n0, Coords *n3,
     double x, double y, double z,
     double kx, double ky, double kz,
     off_struct data )
* ACTION: computes intersection of an xray trajectory with an object.
* INPUT:  x,y,z and kx,ky,kz, are spatial coordinates and wavevector of the x-ray
*         respectively. data points to the OFF data structure.
* RETURN: the number of polyhedra the trajectory intersects
*         l0 and l3 are the smallest incoming and outgoing intersection lengths
*         n0 and n3 are the corresponding normal vectors to the surface
*******************************************************************************/
int off_x_intersect(double *l0,double *l3,
     Coords *n0, Coords *n3,
     double x,  double y,  double z,
     double kx, double ky, double kz,
     off_struct data )
{
  /*This function simply reformats and calls off_intersect (as for neutrons)
   *by normalizing the wavevector - this will yield the intersection lengths
   *in m*/
  double jx,jy,jz,invk;
  int n;
  invk=1/sqrt(scalar_prod(kx,ky,kz,kx,ky,kz));
  jx=kx*invk;jy=ky*invk;jz=kz*invk;
  n=off_intersect(l0,l3,n0,n3,x,y,z,jx,jy,jz,0.0,0.0,0.0,data);
  return n;
}


/*******************************************************************************
* void off_display(off_struct data)
* ACTION: display up to N_VERTEX_DISPLAYED polygons from the object
*******************************************************************************/
void off_display(off_struct data)
{
  unsigned int i;
  double ratio=(double)(N_VERTEX_DISPLAYED)/(double)data.faceSize;
  unsigned int pixel=0;
  for (i=0; i<data.faceSize-1; i++) {
    int j;
    int nbVertex = data.faceArray[i];
    double x0,y0,z0;
    x0 = data.vtxArray[data.faceArray[i+1]].x;
    y0 = data.vtxArray[data.faceArray[i+1]].y;
    z0 = data.vtxArray[data.faceArray[i+1]].z;
    double x1=x0,y1=y0,z1=z0;
    double cmx=0,cmy=0,cmz=0;

    int drawthis = rand01() < ratio;
    // First pass, calculate center of mass location...
    for (j=1; j<=nbVertex; j++) {
      cmx = cmx+data.vtxArray[data.faceArray[i+j]].x;
      cmy = cmy+data.vtxArray[data.faceArray[i+j]].y;
      cmz = cmz+data.vtxArray[data.faceArray[i+j]].z;
    }
    cmx /= nbVertex;
    cmy /= nbVertex;
    cmz /= nbVertex;

    char pixelinfo[1024];
    sprintf(pixelinfo, "%li,%li,%li,%i,%g,%g,%g,%g,%g,%g", data.mantidoffset+pixel, data.mantidoffset, data.mantidoffset+data.polySize-1, nbVertex, cmx, cmy, cmz, x1-cmx, y1-cmy, z1-cmz);
    for (j=2; j<=nbVertex; j++) {
      double x2,y2,z2;
      x2 = data.vtxArray[data.faceArray[i+j]].x;
      y2 = data.vtxArray[data.faceArray[i+j]].y;
      z2 = data.vtxArray[data.faceArray[i+j]].z;
      sprintf(pixelinfo, "%s,%g,%g,%g", pixelinfo, x2-cmx, y2-cmy, z2-cmz);
      if (ratio > 1 || drawthis) {
	mcdis_line(x1,y1,z1,x2,y2,z2);
      }
      x1 = x2; y1 = y2; z1 = z2;
    }
    if (ratio > 1 || drawthis) {
	mcdis_line(x1,y1,z1,x0,y0,z0);
      }
    if (data.mantidflag) {
      printf("MANTID_PIXEL: %s\n", pixelinfo);
      pixel++;
    }
    i += nbVertex;
  }
} /* off_display */

/* end of interoff-lib.c */




/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component origin=Progress_bar() [1] DECLARE */
/* Parameter definition for component type 'Progress_bar' */
struct _struct_Progress_bar_parameters {
  /* Component type 'Progress_bar' setting parameters */
  char profile[16384];
  MCNUM percent;
  MCNUM flag_save;
  MCNUM minutes;
  /* Component type 'Progress_bar' private parameters */
  double  IntermediateCnts;
  time_t  StartTime;
  time_t  EndTime;
  time_t  CurrentTime;
}; /* _struct_Progress_bar_parameters */
typedef struct _struct_Progress_bar_parameters _class_Progress_bar_parameters;

/* Parameters for component type 'Progress_bar' */
struct _struct_Progress_bar {
  char     _name[256]; /* e.g. origin */
  char     _type[256]; /* Progress_bar */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Progress_bar_parameters _parameters;
};
typedef struct _struct_Progress_bar _class_Progress_bar;
_class_Progress_bar _origin_var;
#pragma acc declare create ( _origin_var )

/* component source=Source_div() [2] DECLARE */
/* Parameter definition for component type 'Source_div' */
struct _struct_Source_div_parameters {
  /* Component type 'Source_div' setting parameters */
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM focus_aw;
  MCNUM focus_ah;
  MCNUM E0;
  MCNUM dE;
  MCNUM lambda0;
  MCNUM dlambda;
  MCNUM gauss;
  MCNUM flux;
  /* Component type 'Source_div' private parameters */
  double  sigmah;
  double  sigmav;
  double  p_init;
  double  dist;
  double  focus_xw;
  double  focus_yh;
}; /* _struct_Source_div_parameters */
typedef struct _struct_Source_div_parameters _class_Source_div_parameters;

/* Parameters for component type 'Source_div' */
struct _struct_Source_div {
  char     _name[256]; /* e.g. source */
  char     _type[256]; /* Source_div */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Source_div_parameters _parameters;
};
typedef struct _struct_Source_div _class_Source_div;
_class_Source_div _source_var;
#pragma acc declare create ( _source_var )

/* component psd_test=PSD_monitor() [3] DECLARE */
/* Parameter definition for component type 'PSD_monitor' */
struct _struct_PSD_monitor_parameters {
  /* Component type 'PSD_monitor' setting parameters */
  MCNUM nx;
  MCNUM ny;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM restore_neutron;
  int nowritefile;
  /* Component type 'PSD_monitor' private parameters */
  DArray2d  PSD_N;
  DArray2d  PSD_p;
  DArray2d  PSD_p2;
}; /* _struct_PSD_monitor_parameters */
typedef struct _struct_PSD_monitor_parameters _class_PSD_monitor_parameters;

/* Parameters for component type 'PSD_monitor' */
struct _struct_PSD_monitor {
  char     _name[256]; /* e.g. psd_test */
  char     _type[256]; /* PSD_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_PSD_monitor_parameters _parameters;
};
typedef struct _struct_PSD_monitor _class_PSD_monitor;
_class_PSD_monitor _psd_test_var;
#pragma acc declare create ( _psd_test_var )

/* component e_monitor=E_monitor() [4] DECLARE */
/* Parameter definition for component type 'E_monitor' */
struct _struct_E_monitor_parameters {
  /* Component type 'E_monitor' setting parameters */
  MCNUM nE;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  int nowritefile;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM Emin;
  MCNUM Emax;
  MCNUM restore_neutron;
  /* Component type 'E_monitor' private parameters */
  DArray1d  E_N;
  DArray1d  E_p;
  DArray1d  E_p2;
  double  S_p;
  double  S_pE;
  double  S_pE2;
}; /* _struct_E_monitor_parameters */
typedef struct _struct_E_monitor_parameters _class_E_monitor_parameters;

/* Parameters for component type 'E_monitor' */
struct _struct_E_monitor {
  char     _name[256]; /* e.g. e_monitor */
  char     _type[256]; /* E_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_E_monitor_parameters _parameters;
};
typedef struct _struct_E_monitor _class_E_monitor;
_class_E_monitor _e_monitor_var;
#pragma acc declare create ( _e_monitor_var )

/* component Analyzer_0=Monochromator_bent() [5] DECLARE */
/* Parameter definition for component type 'Monochromator_bent' */
struct _struct_Monochromator_bent_parameters {
  /* Component type 'Monochromator_bent' setting parameters */
  MCNUM zwidth;
  MCNUM yheight;
  MCNUM xthickness;
  MCNUM radius_x;
  char plane_of_reflection[16384];
  MCNUM angle_to_cut_horizontal;
  MCNUM angle_to_cut_vertical;
  MCNUM mosaicity;
  MCNUM domainthickness;
  MCNUM temperature;
  int verbose;
  /* Component type 'Monochromator_bent' private parameters */
  struct neutron_values  neutron;
  struct monochromator_values  monochromator;
  int  non_scattered;
  int  scattered;
  int  non_hit;
}; /* _struct_Monochromator_bent_parameters */
typedef struct _struct_Monochromator_bent_parameters _class_Monochromator_bent_parameters;

/* Parameters for component type 'Monochromator_bent' */
struct _struct_Monochromator_bent {
  char     _name[256]; /* e.g. Analyzer_0 */
  char     _type[256]; /* Monochromator_bent */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Monochromator_bent_parameters _parameters;
};
typedef struct _struct_Monochromator_bent _class_Monochromator_bent;
_class_Monochromator_bent _Analyzer_0_var;
#pragma acc declare create ( _Analyzer_0_var )

_class_Monochromator_bent _Analyzer_1_var;
#pragma acc declare create ( _Analyzer_1_var )

_class_Monochromator_bent _Analyzer_2_var;
#pragma acc declare create ( _Analyzer_2_var )

_class_Monochromator_bent _Analyzer_3_var;
#pragma acc declare create ( _Analyzer_3_var )

_class_Monochromator_bent _Analyzer_4_var;
#pragma acc declare create ( _Analyzer_4_var )

_class_Monochromator_bent _Analyzer_5_var;
#pragma acc declare create ( _Analyzer_5_var )

_class_Monochromator_bent _Analyzer_6_var;
#pragma acc declare create ( _Analyzer_6_var )

_class_Monochromator_bent _Analyzer_7_var;
#pragma acc declare create ( _Analyzer_7_var )

_class_Monochromator_bent _Analyzer_8_var;
#pragma acc declare create ( _Analyzer_8_var )

_class_Monochromator_bent _Analyzer_9_var;
#pragma acc declare create ( _Analyzer_9_var )

_class_Monochromator_bent _Analyzer_10_var;
#pragma acc declare create ( _Analyzer_10_var )

_class_Monochromator_bent _Analyzer_11_var;
#pragma acc declare create ( _Analyzer_11_var )

_class_Monochromator_bent _Analyzer_12_var;
#pragma acc declare create ( _Analyzer_12_var )

_class_Monochromator_bent _Analyzer_13_var;
#pragma acc declare create ( _Analyzer_13_var )

_class_Monochromator_bent _Analyzer_14_var;
#pragma acc declare create ( _Analyzer_14_var )

_class_Monochromator_bent _Analyzer_15_var;
#pragma acc declare create ( _Analyzer_15_var )

_class_Monochromator_bent _Analyzer_16_var;
#pragma acc declare create ( _Analyzer_16_var )

_class_Monochromator_bent _Analyzer_17_var;
#pragma acc declare create ( _Analyzer_17_var )

_class_Monochromator_bent _Analyzer_18_var;
#pragma acc declare create ( _Analyzer_18_var )

_class_Monochromator_bent _Analyzer_19_var;
#pragma acc declare create ( _Analyzer_19_var )

_class_Monochromator_bent _Analyzer_20_var;
#pragma acc declare create ( _Analyzer_20_var )

_class_Monochromator_bent _Analyzer_21_var;
#pragma acc declare create ( _Analyzer_21_var )

_class_Monochromator_bent _Analyzer_22_var;
#pragma acc declare create ( _Analyzer_22_var )

_class_Monochromator_bent _Analyzer_23_var;
#pragma acc declare create ( _Analyzer_23_var )

_class_Monochromator_bent _Analyzer_24_var;
#pragma acc declare create ( _Analyzer_24_var )

_class_Monochromator_bent _Analyzer_25_var;
#pragma acc declare create ( _Analyzer_25_var )

_class_Monochromator_bent _Analyzer_26_var;
#pragma acc declare create ( _Analyzer_26_var )

_class_Monochromator_bent _Analyzer_27_var;
#pragma acc declare create ( _Analyzer_27_var )

_class_Monochromator_bent _Analyzer_28_var;
#pragma acc declare create ( _Analyzer_28_var )

_class_Monochromator_bent _Analyzer_29_var;
#pragma acc declare create ( _Analyzer_29_var )

_class_Monochromator_bent _Analyzer_30_var;
#pragma acc declare create ( _Analyzer_30_var )

_class_Monochromator_bent _Analyzer_31_var;
#pragma acc declare create ( _Analyzer_31_var )

_class_Monochromator_bent _Analyzer_32_var;
#pragma acc declare create ( _Analyzer_32_var )

_class_Monochromator_bent _Analyzer_33_var;
#pragma acc declare create ( _Analyzer_33_var )

_class_Monochromator_bent _Analyzer_34_var;
#pragma acc declare create ( _Analyzer_34_var )

_class_Monochromator_bent _Analyzer_35_var;
#pragma acc declare create ( _Analyzer_35_var )

_class_Monochromator_bent _Analyzer_36_var;
#pragma acc declare create ( _Analyzer_36_var )

_class_Monochromator_bent _Analyzer_37_var;
#pragma acc declare create ( _Analyzer_37_var )

_class_Monochromator_bent _Analyzer_38_var;
#pragma acc declare create ( _Analyzer_38_var )

_class_Monochromator_bent _Analyzer_39_var;
#pragma acc declare create ( _Analyzer_39_var )

_class_Monochromator_bent _Analyzer_40_var;
#pragma acc declare create ( _Analyzer_40_var )

_class_Monochromator_bent _Analyzer_41_var;
#pragma acc declare create ( _Analyzer_41_var )

_class_Monochromator_bent _Analyzer_42_var;
#pragma acc declare create ( _Analyzer_42_var )

_class_Monochromator_bent _Analyzer_43_var;
#pragma acc declare create ( _Analyzer_43_var )

_class_Monochromator_bent _Analyzer_44_var;
#pragma acc declare create ( _Analyzer_44_var )

_class_Monochromator_bent _Analyzer_45_var;
#pragma acc declare create ( _Analyzer_45_var )

/* component rotator=Arm() [51] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. rotator */
  char     _type[256]; /* Arm */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Arm_parameters _parameters;
};
typedef struct _struct_Arm _class_Arm;
_class_Arm _rotator_var;
#pragma acc declare create ( _rotator_var )

_class_Arm _detector_pos_var;
#pragma acc declare create ( _detector_pos_var )

_class_Arm _detector_var;
#pragma acc declare create ( _detector_var )

_class_PSD_monitor _psd_monitor_end_var;
#pragma acc declare create ( _psd_monitor_end_var )

/* component E_PSD_mon_end=Monitor_nD() [55] DECLARE */
/* Parameter definition for component type 'Monitor_nD' */
struct _struct_Monitor_nD_parameters {
  /* Component type 'Monitor_nD' setting parameters */
  char user1[16384];
  char user2[16384];
  char user3[16384];
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM zdepth;
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM zmin;
  MCNUM zmax;
  MCNUM bins;
  MCNUM min;
  MCNUM max;
  MCNUM restore_neutron;
  MCNUM radius;
  char options[16384];
  char filename[16384];
  char geometry[16384];
  int nowritefile;
  char username1[16384];
  char username2[16384];
  char username3[16384];
  /* Component type 'Monitor_nD' private parameters */
  MonitornD_Defines_type  DEFS;
  MonitornD_Variables_type  Vars;
  MCDETECTOR  detector;
  off_struct  offdata;
}; /* _struct_Monitor_nD_parameters */
typedef struct _struct_Monitor_nD_parameters _class_Monitor_nD_parameters;

/* Parameters for component type 'Monitor_nD' */
struct _struct_Monitor_nD {
  char     _name[256]; /* e.g. E_PSD_mon_end */
  char     _type[256]; /* Monitor_nD */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Monitor_nD_parameters _parameters;
};
typedef struct _struct_Monitor_nD _class_Monitor_nD;
_class_Monitor_nD _E_PSD_mon_end_var;
#pragma acc declare create ( _E_PSD_mon_end_var )

int mcNUMCOMP = 55;

/* User declarations from instrument definition. Can define functions. */
int scat=0;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'template_simple' and components DECLARE */

/* *****************************************************************************
* instrument 'template_simple' and components INITIALISE
***************************************************************************** */

double index_getdistance(int first_index, int second_index)
/* Calculate the distance two components from their indexes*/
{
  return coords_len(coords_sub(POS_A_COMP_INDEX(first_index), POS_A_COMP_INDEX(second_index)));
}

double getdistance(char* first_component, char* second_component)
/* Calculate the distance between two named components */
{
  int first_index = _getcomp_index(first_component);
  int second_index = _getcomp_index(second_component);
  return index_getdistance(first_index, second_index);
}

double checked_setpos_getdistance(int current_index, char* first_component, char* second_component)
/* Calculate the distance between two named components at *_setpos() time, with component index checking */
{
  int first_index = _getcomp_index(first_component);
  int second_index = _getcomp_index(second_component);
  if (first_index >= current_index || second_index >= current_index) {
    printf("setpos_getdistance can only be used with the names of components before the current one!\n");
    return 0;
  }
  return index_getdistance(first_index, second_index);
}
#define setpos_getdistance(first, second) checked_setpos_getdistance(current_setpos_index, first, second)

/* component origin=Progress_bar() SETTING, POSITION/ROTATION */
int _origin_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_origin_setpos] component origin=Progress_bar() SETTING [C:\\mcstas-3.4\\lib\\misc\\Progress_bar.comp:57]");
  stracpy(_origin_var._name, "origin", 16384);
  stracpy(_origin_var._type, "Progress_bar", 16384);
  _origin_var._index=1;
  int current_setpos_index = 1;
  if("NULL" && strlen("NULL"))
    stracpy(_origin_var._parameters.profile, "NULL" ? "NULL" : "", 16384);
  else 
  _origin_var._parameters.profile[0]='\0';
  _origin_var._parameters.percent = 10;
  _origin_var._parameters.flag_save = 0;
  _origin_var._parameters.minutes = 0;


  /* component origin=Progress_bar() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(_origin_var._rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_origin_var._rotation_relative, _origin_var._rotation_absolute);
    _origin_var._rotation_is_identity =  rot_test_identity(_origin_var._rotation_relative);
    _origin_var._position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_origin_var._position_absolute);
    _origin_var._position_relative = rot_apply(_origin_var._rotation_absolute, tc1);
  } /* origin=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("origin", _origin_var._position_absolute, _origin_var._rotation_absolute);
  instrument->_position_absolute[1] = _origin_var._position_absolute;
  instrument->_position_relative[1] = _origin_var._position_relative;
    _origin_var._position_relative_is_zero =  coords_test_zero(_origin_var._position_relative);
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _origin_setpos */

/* component source=Source_div() SETTING, POSITION/ROTATION */
int _source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_source_setpos] component source=Source_div() SETTING [C:\\mcstas-3.4\\lib\\sources\\Source_div.comp:77]");
  stracpy(_source_var._name, "source", 16384);
  stracpy(_source_var._type, "Source_div", 16384);
  _source_var._index=2;
  int current_setpos_index = 2;
  _source_var._parameters.xwidth = _instrument_var._parameters.sample_y;
  _source_var._parameters.yheight = _instrument_var._parameters.sample_x;
  _source_var._parameters.focus_aw = 4;
  _source_var._parameters.focus_ah = 2.5;
  _source_var._parameters.E0 = _instrument_var._parameters.E_m;
  _source_var._parameters.dE = _instrument_var._parameters.d_E;
  _source_var._parameters.lambda0 = 0.0;
  _source_var._parameters.dlambda = 0.0;
  _source_var._parameters.gauss = 0;
  _source_var._parameters.flux = 1;


  /* component source=Source_div() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _origin_var._rotation_absolute, _source_var._rotation_absolute);
    rot_transpose(_origin_var._rotation_absolute, tr1);
    rot_mul(_source_var._rotation_absolute, tr1, _source_var._rotation_relative);
    _source_var._rotation_is_identity =  rot_test_identity(_source_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_origin_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _source_var._position_absolute = coords_add(_origin_var._position_absolute, tc2);
    tc1 = coords_sub(_origin_var._position_absolute, _source_var._position_absolute);
    _source_var._position_relative = rot_apply(_source_var._rotation_absolute, tc1);
  } /* source=Source_div() AT ROTATED */
  DEBUG_COMPONENT("source", _source_var._position_absolute, _source_var._rotation_absolute);
  instrument->_position_absolute[2] = _source_var._position_absolute;
  instrument->_position_relative[2] = _source_var._position_relative;
    _source_var._position_relative_is_zero =  coords_test_zero(_source_var._position_relative);
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _source_setpos */

/* component psd_test=PSD_monitor() SETTING, POSITION/ROTATION */
int _psd_test_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psd_test_setpos] component psd_test=PSD_monitor() SETTING [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:62]");
  stracpy(_psd_test_var._name, "psd_test", 16384);
  stracpy(_psd_test_var._type, "PSD_monitor", 16384);
  _psd_test_var._index=3;
  int current_setpos_index = 3;
  _psd_test_var._parameters.nx = 90;
  _psd_test_var._parameters.ny = 90;
  if("psd_test" && strlen("psd_test"))
    stracpy(_psd_test_var._parameters.filename, "psd_test" ? "psd_test" : "", 16384);
  else 
  _psd_test_var._parameters.filename[0]='\0';
  _psd_test_var._parameters.xmin = -0.1;
  _psd_test_var._parameters.xmax = 0.1;
  _psd_test_var._parameters.ymin = -0.1;
  _psd_test_var._parameters.ymax = 0.1;
  _psd_test_var._parameters.xwidth = 0;
  _psd_test_var._parameters.yheight = 0;
  _psd_test_var._parameters.restore_neutron = 1;
  _psd_test_var._parameters.nowritefile = 0;


  /* component psd_test=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _psd_test_var._rotation_absolute);
    rot_transpose(_source_var._rotation_absolute, tr1);
    rot_mul(_psd_test_var._rotation_absolute, tr1, _psd_test_var._rotation_relative);
    _psd_test_var._rotation_is_identity =  rot_test_identity(_psd_test_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psd_test_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_source_var._position_absolute, _psd_test_var._position_absolute);
    _psd_test_var._position_relative = rot_apply(_psd_test_var._rotation_absolute, tc1);
  } /* psd_test=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psd_test", _psd_test_var._position_absolute, _psd_test_var._rotation_absolute);
  instrument->_position_absolute[3] = _psd_test_var._position_absolute;
  instrument->_position_relative[3] = _psd_test_var._position_relative;
    _psd_test_var._position_relative_is_zero =  coords_test_zero(_psd_test_var._position_relative);
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _psd_test_setpos */

/* component e_monitor=E_monitor() SETTING, POSITION/ROTATION */
int _e_monitor_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_e_monitor_setpos] component e_monitor=E_monitor() SETTING [C:\\mcstas-3.4\\lib\\monitors\\E_monitor.comp:69]");
  stracpy(_e_monitor_var._name, "e_monitor", 16384);
  stracpy(_e_monitor_var._type, "E_monitor", 16384);
  _e_monitor_var._index=4;
  int current_setpos_index = 4;
  _e_monitor_var._parameters.nE = 501;
  if("E_start" && strlen("E_start"))
    stracpy(_e_monitor_var._parameters.filename, "E_start" ? "E_start" : "", 16384);
  else 
  _e_monitor_var._parameters.filename[0]='\0';
  _e_monitor_var._parameters.xmin = -0.05;
  _e_monitor_var._parameters.xmax = 0.05;
  _e_monitor_var._parameters.ymin = -0.05;
  _e_monitor_var._parameters.ymax = 0.05;
  _e_monitor_var._parameters.nowritefile = 0;
  _e_monitor_var._parameters.xwidth = 1;
  _e_monitor_var._parameters.yheight = 1;
  _e_monitor_var._parameters.Emin = _instrument_var._parameters.E_m - _instrument_var._parameters.d_E * 1.1;
  _e_monitor_var._parameters.Emax = _instrument_var._parameters.E_m + _instrument_var._parameters.d_E * 1.1;
  _e_monitor_var._parameters.restore_neutron = 1;


  /* component e_monitor=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _e_monitor_var._rotation_absolute);
    rot_transpose(_psd_test_var._rotation_absolute, tr1);
    rot_mul(_e_monitor_var._rotation_absolute, tr1, _e_monitor_var._rotation_relative);
    _e_monitor_var._rotation_is_identity =  rot_test_identity(_e_monitor_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _e_monitor_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_psd_test_var._position_absolute, _e_monitor_var._position_absolute);
    _e_monitor_var._position_relative = rot_apply(_e_monitor_var._rotation_absolute, tc1);
  } /* e_monitor=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("e_monitor", _e_monitor_var._position_absolute, _e_monitor_var._rotation_absolute);
  instrument->_position_absolute[4] = _e_monitor_var._position_absolute;
  instrument->_position_relative[4] = _e_monitor_var._position_relative;
    _e_monitor_var._position_relative_is_zero =  coords_test_zero(_e_monitor_var._position_relative);
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _e_monitor_setpos */

/* component Analyzer_0=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_0_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_0_setpos] component Analyzer_0=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_0_var._name, "Analyzer_0", 16384);
  stracpy(_Analyzer_0_var._type, "Monochromator_bent", 16384);
  _Analyzer_0_var._index=5;
  int current_setpos_index = 5;
  _Analyzer_0_var._parameters.zwidth = 0.045186153053214745;
  _Analyzer_0_var._parameters.yheight = 0.025256242964165475;
  _Analyzer_0_var._parameters.xthickness = 0.0005;
  _Analyzer_0_var._parameters.radius_x = 2.442094641246816;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_0_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_0_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_0_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_0_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_0_var._parameters.mosaicity = 20;
  _Analyzer_0_var._parameters.domainthickness = 10;
  _Analyzer_0_var._parameters.temperature = 0;
  _Analyzer_0_var._parameters.verbose = 1;


  /* component Analyzer_0=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (34.50078023652336)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_0_var._rotation_absolute);
    rot_transpose(_e_monitor_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_0_var._rotation_absolute, tr1, _Analyzer_0_var._rotation_relative);
    _Analyzer_0_var._rotation_is_identity =  rot_test_identity(_Analyzer_0_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5800884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_0_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_e_monitor_var._position_absolute, _Analyzer_0_var._position_absolute);
    _Analyzer_0_var._position_relative = rot_apply(_Analyzer_0_var._rotation_absolute, tc1);
  } /* Analyzer_0=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_0", _Analyzer_0_var._position_absolute, _Analyzer_0_var._rotation_absolute);
  instrument->_position_absolute[5] = _Analyzer_0_var._position_absolute;
  instrument->_position_relative[5] = _Analyzer_0_var._position_relative;
    _Analyzer_0_var._position_relative_is_zero =  coords_test_zero(_Analyzer_0_var._position_relative);
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _Analyzer_0_setpos */

/* component Analyzer_1=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_1_setpos] component Analyzer_1=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_1_var._name, "Analyzer_1", 16384);
  stracpy(_Analyzer_1_var._type, "Monochromator_bent", 16384);
  _Analyzer_1_var._index=6;
  int current_setpos_index = 6;
  _Analyzer_1_var._parameters.zwidth = 0.04656520220403678;
  _Analyzer_1_var._parameters.yheight = 0.025634711938996;
  _Analyzer_1_var._parameters.xthickness = 0.0005;
  _Analyzer_1_var._parameters.radius_x = 2.480458519804657;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_1_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_1_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_1_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_1_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_1_var._parameters.mosaicity = 20;
  _Analyzer_1_var._parameters.domainthickness = 10;
  _Analyzer_1_var._parameters.temperature = 0;
  _Analyzer_1_var._parameters.verbose = 1;


  /* component Analyzer_1=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (34.739608321525225)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_1_var._rotation_absolute);
    rot_transpose(_Analyzer_0_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_1_var._rotation_absolute, tr1, _Analyzer_1_var._rotation_relative);
    _Analyzer_1_var._rotation_is_identity =  rot_test_identity(_Analyzer_1_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5900884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_1_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_0_var._position_absolute, _Analyzer_1_var._position_absolute);
    _Analyzer_1_var._position_relative = rot_apply(_Analyzer_1_var._rotation_absolute, tc1);
  } /* Analyzer_1=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_1", _Analyzer_1_var._position_absolute, _Analyzer_1_var._rotation_absolute);
  instrument->_position_absolute[6] = _Analyzer_1_var._position_absolute;
  instrument->_position_relative[6] = _Analyzer_1_var._position_relative;
    _Analyzer_1_var._position_relative_is_zero =  coords_test_zero(_Analyzer_1_var._position_relative);
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _Analyzer_1_setpos */

/* component Analyzer_2=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_2_setpos] component Analyzer_2=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_2_var._name, "Analyzer_2", 16384);
  stracpy(_Analyzer_2_var._type, "Monochromator_bent", 16384);
  _Analyzer_2_var._index=7;
  int current_setpos_index = 7;
  _Analyzer_2_var._parameters.zwidth = 0.05153642294323662;
  _Analyzer_2_var._parameters.yheight = 0.025941422423863977;
  _Analyzer_2_var._parameters.xthickness = 0.0005;
  _Analyzer_2_var._parameters.radius_x = 2.504311157624173;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_2_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_2_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_2_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_2_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_2_var._parameters.mosaicity = 20;
  _Analyzer_2_var._parameters.domainthickness = 10;
  _Analyzer_2_var._parameters.temperature = 0;
  _Analyzer_2_var._parameters.verbose = 1;


  /* component Analyzer_2=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (34.97302867718735)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_2_var._rotation_absolute);
    rot_transpose(_Analyzer_1_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_2_var._rotation_absolute, tr1, _Analyzer_2_var._rotation_relative);
    _Analyzer_2_var._rotation_is_identity =  rot_test_identity(_Analyzer_2_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6000884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_2_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_1_var._position_absolute, _Analyzer_2_var._position_absolute);
    _Analyzer_2_var._position_relative = rot_apply(_Analyzer_2_var._rotation_absolute, tc1);
  } /* Analyzer_2=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_2", _Analyzer_2_var._position_absolute, _Analyzer_2_var._rotation_absolute);
  instrument->_position_absolute[7] = _Analyzer_2_var._position_absolute;
  instrument->_position_relative[7] = _Analyzer_2_var._position_relative;
    _Analyzer_2_var._position_relative_is_zero =  coords_test_zero(_Analyzer_2_var._position_relative);
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _Analyzer_2_setpos */

/* component Analyzer_3=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_3_setpos] component Analyzer_3=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_3_var._name, "Analyzer_3", 16384);
  stracpy(_Analyzer_3_var._type, "Monochromator_bent", 16384);
  _Analyzer_3_var._index=8;
  int current_setpos_index = 8;
  _Analyzer_3_var._parameters.zwidth = 0.0548909655519918;
  _Analyzer_3_var._parameters.yheight = 0.02624997411290457;
  _Analyzer_3_var._parameters.xthickness = 0.0005;
  _Analyzer_3_var._parameters.radius_x = 2.5092331783495223;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_3_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_3_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_3_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_3_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_3_var._parameters.mosaicity = 20;
  _Analyzer_3_var._parameters.domainthickness = 10;
  _Analyzer_3_var._parameters.temperature = 0;
  _Analyzer_3_var._parameters.verbose = 1;


  /* component Analyzer_3=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (35.20433658263342)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_3_var._rotation_absolute);
    rot_transpose(_Analyzer_2_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_3_var._rotation_absolute, tr1, _Analyzer_3_var._rotation_relative);
    _Analyzer_3_var._rotation_is_identity =  rot_test_identity(_Analyzer_3_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6100884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_3_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_2_var._position_absolute, _Analyzer_3_var._position_absolute);
    _Analyzer_3_var._position_relative = rot_apply(_Analyzer_3_var._rotation_absolute, tc1);
  } /* Analyzer_3=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_3", _Analyzer_3_var._position_absolute, _Analyzer_3_var._rotation_absolute);
  instrument->_position_absolute[8] = _Analyzer_3_var._position_absolute;
  instrument->_position_relative[8] = _Analyzer_3_var._position_relative;
    _Analyzer_3_var._position_relative_is_zero =  coords_test_zero(_Analyzer_3_var._position_relative);
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _Analyzer_3_setpos */

/* component Analyzer_4=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_4_setpos] component Analyzer_4=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_4_var._name, "Analyzer_4", 16384);
  stracpy(_Analyzer_4_var._type, "Monochromator_bent", 16384);
  _Analyzer_4_var._index=9;
  int current_setpos_index = 9;
  _Analyzer_4_var._parameters.zwidth = 0.05740910355091564;
  _Analyzer_4_var._parameters.yheight = 0.026589164467988584;
  _Analyzer_4_var._parameters.xthickness = 0.0005;
  _Analyzer_4_var._parameters.radius_x = 2.527981711977253;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_4_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_4_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_4_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_4_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_4_var._parameters.mosaicity = 20;
  _Analyzer_4_var._parameters.domainthickness = 10;
  _Analyzer_4_var._parameters.temperature = 0;
  _Analyzer_4_var._parameters.verbose = 1;


  /* component Analyzer_4=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (35.438817405922826)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_4_var._rotation_absolute);
    rot_transpose(_Analyzer_3_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_4_var._rotation_absolute, tr1, _Analyzer_4_var._rotation_relative);
    _Analyzer_4_var._rotation_is_identity =  rot_test_identity(_Analyzer_4_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6200884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_4_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_3_var._position_absolute, _Analyzer_4_var._position_absolute);
    _Analyzer_4_var._position_relative = rot_apply(_Analyzer_4_var._rotation_absolute, tc1);
  } /* Analyzer_4=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_4", _Analyzer_4_var._position_absolute, _Analyzer_4_var._rotation_absolute);
  instrument->_position_absolute[9] = _Analyzer_4_var._position_absolute;
  instrument->_position_relative[9] = _Analyzer_4_var._position_relative;
    _Analyzer_4_var._position_relative_is_zero =  coords_test_zero(_Analyzer_4_var._position_relative);
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _Analyzer_4_setpos */

/* component Analyzer_5=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_5_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_5_setpos] component Analyzer_5=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_5_var._name, "Analyzer_5", 16384);
  stracpy(_Analyzer_5_var._type, "Monochromator_bent", 16384);
  _Analyzer_5_var._index=10;
  int current_setpos_index = 10;
  _Analyzer_5_var._parameters.zwidth = 0.06144350178367934;
  _Analyzer_5_var._parameters.yheight = 0.02694414833540672;
  _Analyzer_5_var._parameters.xthickness = 0.0005;
  _Analyzer_5_var._parameters.radius_x = 2.5742328029536248;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_5_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_5_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_5_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_5_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_5_var._parameters.mosaicity = 20;
  _Analyzer_5_var._parameters.domainthickness = 10;
  _Analyzer_5_var._parameters.temperature = 0;
  _Analyzer_5_var._parameters.verbose = 1;


  /* component Analyzer_5=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (35.67618426233606)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_5_var._rotation_absolute);
    rot_transpose(_Analyzer_4_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_5_var._rotation_absolute, tr1, _Analyzer_5_var._rotation_relative);
    _Analyzer_5_var._rotation_is_identity =  rot_test_identity(_Analyzer_5_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6300884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_5_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_4_var._position_absolute, _Analyzer_5_var._position_absolute);
    _Analyzer_5_var._position_relative = rot_apply(_Analyzer_5_var._rotation_absolute, tc1);
  } /* Analyzer_5=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_5", _Analyzer_5_var._position_absolute, _Analyzer_5_var._rotation_absolute);
  instrument->_position_absolute[10] = _Analyzer_5_var._position_absolute;
  instrument->_position_relative[10] = _Analyzer_5_var._position_relative;
    _Analyzer_5_var._position_relative_is_zero =  coords_test_zero(_Analyzer_5_var._position_relative);
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _Analyzer_5_setpos */

/* component Analyzer_6=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_6_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_6_setpos] component Analyzer_6=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_6_var._name, "Analyzer_6", 16384);
  stracpy(_Analyzer_6_var._type, "Monochromator_bent", 16384);
  _Analyzer_6_var._index=11;
  int current_setpos_index = 11;
  _Analyzer_6_var._parameters.zwidth = 0.06391907703861957;
  _Analyzer_6_var._parameters.yheight = 0.02727198937355823;
  _Analyzer_6_var._parameters.xthickness = 0.0005;
  _Analyzer_6_var._parameters.radius_x = 2.5825867906391307;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_6_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_6_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_6_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_6_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_6_var._parameters.mosaicity = 20;
  _Analyzer_6_var._parameters.domainthickness = 10;
  _Analyzer_6_var._parameters.temperature = 0;
  _Analyzer_6_var._parameters.verbose = 1;


  /* component Analyzer_6=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (35.909330306732805)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_6_var._rotation_absolute);
    rot_transpose(_Analyzer_5_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_6_var._rotation_absolute, tr1, _Analyzer_6_var._rotation_relative);
    _Analyzer_6_var._rotation_is_identity =  rot_test_identity(_Analyzer_6_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6400884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_6_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_5_var._position_absolute, _Analyzer_6_var._position_absolute);
    _Analyzer_6_var._position_relative = rot_apply(_Analyzer_6_var._rotation_absolute, tc1);
  } /* Analyzer_6=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_6", _Analyzer_6_var._position_absolute, _Analyzer_6_var._rotation_absolute);
  instrument->_position_absolute[11] = _Analyzer_6_var._position_absolute;
  instrument->_position_relative[11] = _Analyzer_6_var._position_relative;
    _Analyzer_6_var._position_relative_is_zero =  coords_test_zero(_Analyzer_6_var._position_relative);
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _Analyzer_6_setpos */

/* component Analyzer_7=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_7_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_7_setpos] component Analyzer_7=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_7_var._name, "Analyzer_7", 16384);
  stracpy(_Analyzer_7_var._type, "Monochromator_bent", 16384);
  _Analyzer_7_var._index=12;
  int current_setpos_index = 12;
  _Analyzer_7_var._parameters.zwidth = 0.06715434288061861;
  _Analyzer_7_var._parameters.yheight = 0.027587314734483063;
  _Analyzer_7_var._parameters.xthickness = 0.0005;
  _Analyzer_7_var._parameters.radius_x = 2.5898657643794;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_7_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_7_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_7_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_7_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_7_var._parameters.mosaicity = 20;
  _Analyzer_7_var._parameters.domainthickness = 10;
  _Analyzer_7_var._parameters.temperature = 0;
  _Analyzer_7_var._parameters.verbose = 1;


  /* component Analyzer_7=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (36.14234958070264)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_7_var._rotation_absolute);
    rot_transpose(_Analyzer_6_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_7_var._rotation_absolute, tr1, _Analyzer_7_var._rotation_relative);
    _Analyzer_7_var._rotation_is_identity =  rot_test_identity(_Analyzer_7_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6500884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_7_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_6_var._position_absolute, _Analyzer_7_var._position_absolute);
    _Analyzer_7_var._position_relative = rot_apply(_Analyzer_7_var._rotation_absolute, tc1);
  } /* Analyzer_7=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_7", _Analyzer_7_var._position_absolute, _Analyzer_7_var._rotation_absolute);
  instrument->_position_absolute[12] = _Analyzer_7_var._position_absolute;
  instrument->_position_relative[12] = _Analyzer_7_var._position_relative;
    _Analyzer_7_var._position_relative_is_zero =  coords_test_zero(_Analyzer_7_var._position_relative);
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _Analyzer_7_setpos */

/* component Analyzer_8=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_8_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_8_setpos] component Analyzer_8=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_8_var._name, "Analyzer_8", 16384);
  stracpy(_Analyzer_8_var._type, "Monochromator_bent", 16384);
  _Analyzer_8_var._index=13;
  int current_setpos_index = 13;
  _Analyzer_8_var._parameters.zwidth = 0.07227180425245408;
  _Analyzer_8_var._parameters.yheight = 0.027904412222491095;
  _Analyzer_8_var._parameters.xthickness = 0.0005;
  _Analyzer_8_var._parameters.radius_x = 2.6161644898079635;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_8_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_8_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_8_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_8_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_8_var._parameters.mosaicity = 20;
  _Analyzer_8_var._parameters.domainthickness = 10;
  _Analyzer_8_var._parameters.temperature = 0;
  _Analyzer_8_var._parameters.verbose = 1;


  /* component Analyzer_8=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (36.37646129104227)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_8_var._rotation_absolute);
    rot_transpose(_Analyzer_7_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_8_var._rotation_absolute, tr1, _Analyzer_8_var._rotation_relative);
    _Analyzer_8_var._rotation_is_identity =  rot_test_identity(_Analyzer_8_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6600884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_8_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_7_var._position_absolute, _Analyzer_8_var._position_absolute);
    _Analyzer_8_var._position_relative = rot_apply(_Analyzer_8_var._rotation_absolute, tc1);
  } /* Analyzer_8=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_8", _Analyzer_8_var._position_absolute, _Analyzer_8_var._rotation_absolute);
  instrument->_position_absolute[13] = _Analyzer_8_var._position_absolute;
  instrument->_position_relative[13] = _Analyzer_8_var._position_relative;
    _Analyzer_8_var._position_relative_is_zero =  coords_test_zero(_Analyzer_8_var._position_relative);
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _Analyzer_8_setpos */

/* component Analyzer_9=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_9_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_9_setpos] component Analyzer_9=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_9_var._name, "Analyzer_9", 16384);
  stracpy(_Analyzer_9_var._type, "Monochromator_bent", 16384);
  _Analyzer_9_var._index=14;
  int current_setpos_index = 14;
  _Analyzer_9_var._parameters.zwidth = 0.07504378187172768;
  _Analyzer_9_var._parameters.yheight = 0.028237260558690767;
  _Analyzer_9_var._parameters.xthickness = 0.0005;
  _Analyzer_9_var._parameters.radius_x = 2.6294769277750336;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_9_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_9_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_9_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_9_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_9_var._parameters.mosaicity = 20;
  _Analyzer_9_var._parameters.domainthickness = 10;
  _Analyzer_9_var._parameters.temperature = 0;
  _Analyzer_9_var._parameters.verbose = 1;


  /* component Analyzer_9=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (36.61112439107199)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_9_var._rotation_absolute);
    rot_transpose(_Analyzer_8_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_9_var._rotation_absolute, tr1, _Analyzer_9_var._rotation_relative);
    _Analyzer_9_var._rotation_is_identity =  rot_test_identity(_Analyzer_9_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.6700884582982309);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_9_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_8_var._position_absolute, _Analyzer_9_var._position_absolute);
    _Analyzer_9_var._position_relative = rot_apply(_Analyzer_9_var._rotation_absolute, tc1);
  } /* Analyzer_9=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_9", _Analyzer_9_var._position_absolute, _Analyzer_9_var._rotation_absolute);
  instrument->_position_absolute[14] = _Analyzer_9_var._position_absolute;
  instrument->_position_relative[14] = _Analyzer_9_var._position_relative;
    _Analyzer_9_var._position_relative_is_zero =  coords_test_zero(_Analyzer_9_var._position_relative);
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _Analyzer_9_setpos */

/* component Analyzer_10=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_10_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_10_setpos] component Analyzer_10=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_10_var._name, "Analyzer_10", 16384);
  stracpy(_Analyzer_10_var._type, "Monochromator_bent", 16384);
  _Analyzer_10_var._index=15;
  int current_setpos_index = 15;
  _Analyzer_10_var._parameters.zwidth = 0.07621584471994772;
  _Analyzer_10_var._parameters.yheight = 0.028627252374239767;
  _Analyzer_10_var._parameters.xthickness = 0.0005;
  _Analyzer_10_var._parameters.radius_x = 2.663728902910889;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_10_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_10_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_10_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_10_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_10_var._parameters.mosaicity = 20;
  _Analyzer_10_var._parameters.domainthickness = 10;
  _Analyzer_10_var._parameters.temperature = 0;
  _Analyzer_10_var._parameters.verbose = 1;


  /* component Analyzer_10=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (36.848848954561525)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_10_var._rotation_absolute);
    rot_transpose(_Analyzer_9_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_10_var._rotation_absolute, tr1, _Analyzer_10_var._rotation_relative);
    _Analyzer_10_var._rotation_is_identity =  rot_test_identity(_Analyzer_10_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.680088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_10_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_9_var._position_absolute, _Analyzer_10_var._position_absolute);
    _Analyzer_10_var._position_relative = rot_apply(_Analyzer_10_var._rotation_absolute, tc1);
  } /* Analyzer_10=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_10", _Analyzer_10_var._position_absolute, _Analyzer_10_var._rotation_absolute);
  instrument->_position_absolute[15] = _Analyzer_10_var._position_absolute;
  instrument->_position_relative[15] = _Analyzer_10_var._position_relative;
    _Analyzer_10_var._position_relative_is_zero =  coords_test_zero(_Analyzer_10_var._position_relative);
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _Analyzer_10_setpos */

/* component Analyzer_11=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_11_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_11_setpos] component Analyzer_11=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_11_var._name, "Analyzer_11", 16384);
  stracpy(_Analyzer_11_var._type, "Monochromator_bent", 16384);
  _Analyzer_11_var._index=16;
  int current_setpos_index = 16;
  _Analyzer_11_var._parameters.zwidth = 0.0792736809977705;
  _Analyzer_11_var._parameters.yheight = 0.02900427151033862;
  _Analyzer_11_var._parameters.xthickness = 0.0005;
  _Analyzer_11_var._parameters.radius_x = 2.7033831580401935;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_11_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_11_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_11_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_11_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_11_var._parameters.mosaicity = 20;
  _Analyzer_11_var._parameters.domainthickness = 10;
  _Analyzer_11_var._parameters.temperature = 0;
  _Analyzer_11_var._parameters.verbose = 1;


  /* component Analyzer_11=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (37.0857658084779)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_11_var._rotation_absolute);
    rot_transpose(_Analyzer_10_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_11_var._rotation_absolute, tr1, _Analyzer_11_var._rotation_relative);
    _Analyzer_11_var._rotation_is_identity =  rot_test_identity(_Analyzer_11_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.690088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_11_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_10_var._position_absolute, _Analyzer_11_var._position_absolute);
    _Analyzer_11_var._position_relative = rot_apply(_Analyzer_11_var._rotation_absolute, tc1);
  } /* Analyzer_11=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_11", _Analyzer_11_var._position_absolute, _Analyzer_11_var._rotation_absolute);
  instrument->_position_absolute[16] = _Analyzer_11_var._position_absolute;
  instrument->_position_relative[16] = _Analyzer_11_var._position_relative;
    _Analyzer_11_var._position_relative_is_zero =  coords_test_zero(_Analyzer_11_var._position_relative);
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _Analyzer_11_setpos */

/* component Analyzer_12=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_12_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_12_setpos] component Analyzer_12=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_12_var._name, "Analyzer_12", 16384);
  stracpy(_Analyzer_12_var._type, "Monochromator_bent", 16384);
  _Analyzer_12_var._index=17;
  int current_setpos_index = 17;
  _Analyzer_12_var._parameters.zwidth = 0.07885720278248978;
  _Analyzer_12_var._parameters.yheight = 0.029450552243892078;
  _Analyzer_12_var._parameters.xthickness = 0.0005;
  _Analyzer_12_var._parameters.radius_x = 2.7510739791690146;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_12_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_12_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_12_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_12_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_12_var._parameters.mosaicity = 20;
  _Analyzer_12_var._parameters.domainthickness = 10;
  _Analyzer_12_var._parameters.temperature = 0;
  _Analyzer_12_var._parameters.verbose = 1;


  /* component Analyzer_12=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (37.32525975037311)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_12_var._rotation_absolute);
    rot_transpose(_Analyzer_11_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_12_var._rotation_absolute, tr1, _Analyzer_12_var._rotation_relative);
    _Analyzer_12_var._rotation_is_identity =  rot_test_identity(_Analyzer_12_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.700088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_12_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_11_var._position_absolute, _Analyzer_12_var._position_absolute);
    _Analyzer_12_var._position_relative = rot_apply(_Analyzer_12_var._rotation_absolute, tc1);
  } /* Analyzer_12=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_12", _Analyzer_12_var._position_absolute, _Analyzer_12_var._rotation_absolute);
  instrument->_position_absolute[17] = _Analyzer_12_var._position_absolute;
  instrument->_position_relative[17] = _Analyzer_12_var._position_relative;
    _Analyzer_12_var._position_relative_is_zero =  coords_test_zero(_Analyzer_12_var._position_relative);
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _Analyzer_12_setpos */

/* component Analyzer_13=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_13_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_13_setpos] component Analyzer_13=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_13_var._name, "Analyzer_13", 16384);
  stracpy(_Analyzer_13_var._type, "Monochromator_bent", 16384);
  _Analyzer_13_var._index=18;
  int current_setpos_index = 18;
  _Analyzer_13_var._parameters.zwidth = 0.07993738527728973;
  _Analyzer_13_var._parameters.yheight = 0.029896654869128132;
  _Analyzer_13_var._parameters.xthickness = 0.0005;
  _Analyzer_13_var._parameters.radius_x = 2.803575207125065;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_13_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_13_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_13_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_13_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_13_var._parameters.mosaicity = 20;
  _Analyzer_13_var._parameters.domainthickness = 10;
  _Analyzer_13_var._parameters.temperature = 0;
  _Analyzer_13_var._parameters.verbose = 1;


  /* component Analyzer_13=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (37.564246303351275)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_13_var._rotation_absolute);
    rot_transpose(_Analyzer_12_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_13_var._rotation_absolute, tr1, _Analyzer_13_var._rotation_relative);
    _Analyzer_13_var._rotation_is_identity =  rot_test_identity(_Analyzer_13_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.710088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_13_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_12_var._position_absolute, _Analyzer_13_var._position_absolute);
    _Analyzer_13_var._position_relative = rot_apply(_Analyzer_13_var._rotation_absolute, tc1);
  } /* Analyzer_13=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_13", _Analyzer_13_var._position_absolute, _Analyzer_13_var._rotation_absolute);
  instrument->_position_absolute[18] = _Analyzer_13_var._position_absolute;
  instrument->_position_relative[18] = _Analyzer_13_var._position_relative;
    _Analyzer_13_var._position_relative_is_zero =  coords_test_zero(_Analyzer_13_var._position_relative);
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _Analyzer_13_setpos */

/* component Analyzer_14=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_14_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_14_setpos] component Analyzer_14=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_14_var._name, "Analyzer_14", 16384);
  stracpy(_Analyzer_14_var._type, "Monochromator_bent", 16384);
  _Analyzer_14_var._index=19;
  int current_setpos_index = 19;
  _Analyzer_14_var._parameters.zwidth = 0.08106362046918529;
  _Analyzer_14_var._parameters.yheight = 0.03028902926957882;
  _Analyzer_14_var._parameters.xthickness = 0.0005;
  _Analyzer_14_var._parameters.radius_x = 2.8265116461070354;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_14_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_14_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_14_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_14_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_14_var._parameters.mosaicity = 20;
  _Analyzer_14_var._parameters.domainthickness = 10;
  _Analyzer_14_var._parameters.temperature = 0;
  _Analyzer_14_var._parameters.verbose = 1;


  /* component Analyzer_14=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (37.801201037118304)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_14_var._rotation_absolute);
    rot_transpose(_Analyzer_13_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_14_var._rotation_absolute, tr1, _Analyzer_14_var._rotation_relative);
    _Analyzer_14_var._rotation_is_identity =  rot_test_identity(_Analyzer_14_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.720088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_14_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_13_var._position_absolute, _Analyzer_14_var._position_absolute);
    _Analyzer_14_var._position_relative = rot_apply(_Analyzer_14_var._rotation_absolute, tc1);
  } /* Analyzer_14=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_14", _Analyzer_14_var._position_absolute, _Analyzer_14_var._rotation_absolute);
  instrument->_position_absolute[19] = _Analyzer_14_var._position_absolute;
  instrument->_position_relative[19] = _Analyzer_14_var._position_relative;
    _Analyzer_14_var._position_relative_is_zero =  coords_test_zero(_Analyzer_14_var._position_relative);
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _Analyzer_14_setpos */

/* component Analyzer_15=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_15_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_15_setpos] component Analyzer_15=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_15_var._name, "Analyzer_15", 16384);
  stracpy(_Analyzer_15_var._type, "Monochromator_bent", 16384);
  _Analyzer_15_var._index=20;
  int current_setpos_index = 20;
  _Analyzer_15_var._parameters.zwidth = 0.08064847614150639;
  _Analyzer_15_var._parameters.yheight = 0.030735321475835346;
  _Analyzer_15_var._parameters.xthickness = 0.0005;
  _Analyzer_15_var._parameters.radius_x = 2.8605208792572774;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_15_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_15_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_15_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_15_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_15_var._parameters.mosaicity = 20;
  _Analyzer_15_var._parameters.domainthickness = 10;
  _Analyzer_15_var._parameters.temperature = 0;
  _Analyzer_15_var._parameters.verbose = 1;


  /* component Analyzer_15=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (38.03983478578341)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_15_var._rotation_absolute);
    rot_transpose(_Analyzer_14_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_15_var._rotation_absolute, tr1, _Analyzer_15_var._rotation_relative);
    _Analyzer_15_var._rotation_is_identity =  rot_test_identity(_Analyzer_15_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.730088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_15_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_14_var._position_absolute, _Analyzer_15_var._position_absolute);
    _Analyzer_15_var._position_relative = rot_apply(_Analyzer_15_var._rotation_absolute, tc1);
  } /* Analyzer_15=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_15", _Analyzer_15_var._position_absolute, _Analyzer_15_var._rotation_absolute);
  instrument->_position_absolute[20] = _Analyzer_15_var._position_absolute;
  instrument->_position_relative[20] = _Analyzer_15_var._position_relative;
    _Analyzer_15_var._position_relative_is_zero =  coords_test_zero(_Analyzer_15_var._position_relative);
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _Analyzer_15_setpos */

/* component Analyzer_16=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_16_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_16_setpos] component Analyzer_16=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_16_var._name, "Analyzer_16", 16384);
  stracpy(_Analyzer_16_var._type, "Monochromator_bent", 16384);
  _Analyzer_16_var._index=21;
  int current_setpos_index = 21;
  _Analyzer_16_var._parameters.zwidth = 0.08207365739332309;
  _Analyzer_16_var._parameters.yheight = 0.031181469483854956;
  _Analyzer_16_var._parameters.xthickness = 0.0005;
  _Analyzer_16_var._parameters.radius_x = 2.897429934745591;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_16_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_16_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_16_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_16_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_16_var._parameters.mosaicity = 20;
  _Analyzer_16_var._parameters.domainthickness = 10;
  _Analyzer_16_var._parameters.temperature = 0;
  _Analyzer_16_var._parameters.verbose = 1;


  /* component Analyzer_16=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (38.27840680026829)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_16_var._rotation_absolute);
    rot_transpose(_Analyzer_15_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_16_var._rotation_absolute, tr1, _Analyzer_16_var._rotation_relative);
    _Analyzer_16_var._rotation_is_identity =  rot_test_identity(_Analyzer_16_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.740088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_16_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_15_var._position_absolute, _Analyzer_16_var._position_absolute);
    _Analyzer_16_var._position_relative = rot_apply(_Analyzer_16_var._rotation_absolute, tc1);
  } /* Analyzer_16=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_16", _Analyzer_16_var._position_absolute, _Analyzer_16_var._rotation_absolute);
  instrument->_position_absolute[21] = _Analyzer_16_var._position_absolute;
  instrument->_position_relative[21] = _Analyzer_16_var._position_relative;
    _Analyzer_16_var._position_relative_is_zero =  coords_test_zero(_Analyzer_16_var._position_relative);
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _Analyzer_16_setpos */

/* component Analyzer_17=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_17_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_17_setpos] component Analyzer_17=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_17_var._name, "Analyzer_17", 16384);
  stracpy(_Analyzer_17_var._type, "Monochromator_bent", 16384);
  _Analyzer_17_var._index=22;
  int current_setpos_index = 22;
  _Analyzer_17_var._parameters.zwidth = 0.083172953001956;
  _Analyzer_17_var._parameters.yheight = 0.03157530493878471;
  _Analyzer_17_var._parameters.xthickness = 0.0005;
  _Analyzer_17_var._parameters.radius_x = 2.910458315461202;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_17_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_17_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_17_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_17_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_17_var._parameters.mosaicity = 20;
  _Analyzer_17_var._parameters.domainthickness = 10;
  _Analyzer_17_var._parameters.temperature = 0;
  _Analyzer_17_var._parameters.verbose = 1;


  /* component Analyzer_17=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (38.51593235684239)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_17_var._rotation_absolute);
    rot_transpose(_Analyzer_16_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_17_var._rotation_absolute, tr1, _Analyzer_17_var._rotation_relative);
    _Analyzer_17_var._rotation_is_identity =  rot_test_identity(_Analyzer_17_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.750088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_17_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_16_var._position_absolute, _Analyzer_17_var._position_absolute);
    _Analyzer_17_var._position_relative = rot_apply(_Analyzer_17_var._rotation_absolute, tc1);
  } /* Analyzer_17=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_17", _Analyzer_17_var._position_absolute, _Analyzer_17_var._rotation_absolute);
  instrument->_position_absolute[22] = _Analyzer_17_var._position_absolute;
  instrument->_position_relative[22] = _Analyzer_17_var._position_relative;
    _Analyzer_17_var._position_relative_is_zero =  coords_test_zero(_Analyzer_17_var._position_relative);
  instrument->counter_N[22]  = instrument->counter_P[22] = instrument->counter_P2[22] = 0;
  instrument->counter_AbsorbProp[22]= 0;
  return(0);
} /* _Analyzer_17_setpos */

/* component Analyzer_18=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_18_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_18_setpos] component Analyzer_18=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_18_var._name, "Analyzer_18", 16384);
  stracpy(_Analyzer_18_var._type, "Monochromator_bent", 16384);
  _Analyzer_18_var._index=23;
  int current_setpos_index = 23;
  _Analyzer_18_var._parameters.zwidth = 0.08275558395724007;
  _Analyzer_18_var._parameters.yheight = 0.032021670575887555;
  _Analyzer_18_var._parameters.xthickness = 0.0005;
  _Analyzer_18_var._parameters.radius_x = 2.9299702571461985;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_18_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_18_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_18_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_18_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_18_var._parameters.mosaicity = 20;
  _Analyzer_18_var._parameters.domainthickness = 10;
  _Analyzer_18_var._parameters.temperature = 0;
  _Analyzer_18_var._parameters.verbose = 1;


  /* component Analyzer_18=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (38.75495465968723)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_18_var._rotation_absolute);
    rot_transpose(_Analyzer_17_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_18_var._rotation_absolute, tr1, _Analyzer_18_var._rotation_relative);
    _Analyzer_18_var._rotation_is_identity =  rot_test_identity(_Analyzer_18_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.760088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_18_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_17_var._position_absolute, _Analyzer_18_var._position_absolute);
    _Analyzer_18_var._position_relative = rot_apply(_Analyzer_18_var._rotation_absolute, tc1);
  } /* Analyzer_18=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_18", _Analyzer_18_var._position_absolute, _Analyzer_18_var._rotation_absolute);
  instrument->_position_absolute[23] = _Analyzer_18_var._position_absolute;
  instrument->_position_relative[23] = _Analyzer_18_var._position_relative;
    _Analyzer_18_var._position_relative_is_zero =  coords_test_zero(_Analyzer_18_var._position_relative);
  instrument->counter_N[23]  = instrument->counter_P[23] = instrument->counter_P2[23] = 0;
  instrument->counter_AbsorbProp[23]= 0;
  return(0);
} /* _Analyzer_18_setpos */

/* component Analyzer_19=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_19_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_19_setpos] component Analyzer_19=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_19_var._name, "Analyzer_19", 16384);
  stracpy(_Analyzer_19_var._type, "Monochromator_bent", 16384);
  _Analyzer_19_var._index=24;
  int current_setpos_index = 24;
  _Analyzer_19_var._parameters.zwidth = 0.08378873237209834;
  _Analyzer_19_var._parameters.yheight = 0.0324679162366587;
  _Analyzer_19_var._parameters.xthickness = 0.0005;
  _Analyzer_19_var._parameters.radius_x = 2.946970319388694;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_19_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_19_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_19_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_19_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_19_var._parameters.mosaicity = 20;
  _Analyzer_19_var._parameters.domainthickness = 10;
  _Analyzer_19_var._parameters.temperature = 0;
  _Analyzer_19_var._parameters.verbose = 1;


  /* component Analyzer_19=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (38.994300590374216)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_19_var._rotation_absolute);
    rot_transpose(_Analyzer_18_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_19_var._rotation_absolute, tr1, _Analyzer_19_var._rotation_relative);
    _Analyzer_19_var._rotation_is_identity =  rot_test_identity(_Analyzer_19_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.770088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_19_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_18_var._position_absolute, _Analyzer_19_var._position_absolute);
    _Analyzer_19_var._position_relative = rot_apply(_Analyzer_19_var._rotation_absolute, tc1);
  } /* Analyzer_19=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_19", _Analyzer_19_var._position_absolute, _Analyzer_19_var._rotation_absolute);
  instrument->_position_absolute[24] = _Analyzer_19_var._position_absolute;
  instrument->_position_relative[24] = _Analyzer_19_var._position_relative;
    _Analyzer_19_var._position_relative_is_zero =  coords_test_zero(_Analyzer_19_var._position_relative);
  instrument->counter_N[24]  = instrument->counter_P[24] = instrument->counter_P2[24] = 0;
  instrument->counter_AbsorbProp[24]= 0;
  return(0);
} /* _Analyzer_19_setpos */

/* component Analyzer_20=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_20_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_20_setpos] component Analyzer_20=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_20_var._name, "Analyzer_20", 16384);
  stracpy(_Analyzer_20_var._type, "Monochromator_bent", 16384);
  _Analyzer_20_var._index=25;
  int current_setpos_index = 25;
  _Analyzer_20_var._parameters.zwidth = 0.08523510527692972;
  _Analyzer_20_var._parameters.yheight = 0.03285047399363413;
  _Analyzer_20_var._parameters.xthickness = 0.0005;
  _Analyzer_20_var._parameters.radius_x = 2.947423462062219;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_20_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_20_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_20_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_20_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_20_var._parameters.mosaicity = 20;
  _Analyzer_20_var._parameters.domainthickness = 10;
  _Analyzer_20_var._parameters.temperature = 0;
  _Analyzer_20_var._parameters.verbose = 1;


  /* component Analyzer_20=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (39.2332829929821)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_20_var._rotation_absolute);
    rot_transpose(_Analyzer_19_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_20_var._rotation_absolute, tr1, _Analyzer_20_var._rotation_relative);
    _Analyzer_20_var._rotation_is_identity =  rot_test_identity(_Analyzer_20_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.780088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_20_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_19_var._position_absolute, _Analyzer_20_var._position_absolute);
    _Analyzer_20_var._position_relative = rot_apply(_Analyzer_20_var._rotation_absolute, tc1);
  } /* Analyzer_20=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_20", _Analyzer_20_var._position_absolute, _Analyzer_20_var._rotation_absolute);
  instrument->_position_absolute[25] = _Analyzer_20_var._position_absolute;
  instrument->_position_relative[25] = _Analyzer_20_var._position_relative;
    _Analyzer_20_var._position_relative_is_zero =  coords_test_zero(_Analyzer_20_var._position_relative);
  instrument->counter_N[25]  = instrument->counter_P[25] = instrument->counter_P2[25] = 0;
  instrument->counter_AbsorbProp[25]= 0;
  return(0);
} /* _Analyzer_20_setpos */

/* component Analyzer_21=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_21_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_21_setpos] component Analyzer_21=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_21_var._name, "Analyzer_21", 16384);
  stracpy(_Analyzer_21_var._type, "Monochromator_bent", 16384);
  _Analyzer_21_var._index=26;
  int current_setpos_index = 26;
  _Analyzer_21_var._parameters.zwidth = 0.08481361434789979;
  _Analyzer_21_var._parameters.yheight = 0.03329707148923055;
  _Analyzer_21_var._parameters.xthickness = 0.0005;
  _Analyzer_21_var._parameters.radius_x = 2.9516948607194715;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_21_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_21_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_21_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_21_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_21_var._parameters.mosaicity = 20;
  _Analyzer_21_var._parameters.domainthickness = 10;
  _Analyzer_21_var._parameters.temperature = 0;
  _Analyzer_21_var._parameters.verbose = 1;


  /* component Analyzer_21=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (39.473861408517664)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_21_var._rotation_absolute);
    rot_transpose(_Analyzer_20_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_21_var._rotation_absolute, tr1, _Analyzer_21_var._rotation_relative);
    _Analyzer_21_var._rotation_is_identity =  rot_test_identity(_Analyzer_21_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.790088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_21_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_20_var._position_absolute, _Analyzer_21_var._position_absolute);
    _Analyzer_21_var._position_relative = rot_apply(_Analyzer_21_var._rotation_absolute, tc1);
  } /* Analyzer_21=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_21", _Analyzer_21_var._position_absolute, _Analyzer_21_var._rotation_absolute);
  instrument->_position_absolute[26] = _Analyzer_21_var._position_absolute;
  instrument->_position_relative[26] = _Analyzer_21_var._position_relative;
    _Analyzer_21_var._position_relative_is_zero =  coords_test_zero(_Analyzer_21_var._position_relative);
  instrument->counter_N[26]  = instrument->counter_P[26] = instrument->counter_P2[26] = 0;
  instrument->counter_AbsorbProp[26]= 0;
  return(0);
} /* _Analyzer_21_setpos */

/* component Analyzer_22=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_22_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_22_setpos] component Analyzer_22=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_22_var._name, "Analyzer_22", 16384);
  stracpy(_Analyzer_22_var._type, "Monochromator_bent", 16384);
  _Analyzer_22_var._index=27;
  int current_setpos_index = 27;
  _Analyzer_22_var._parameters.zwidth = 0.08582014507208968;
  _Analyzer_22_var._parameters.yheight = 0.033743566706154735;
  _Analyzer_22_var._parameters.xthickness = 0.0005;
  _Analyzer_22_var._parameters.radius_x = 2.9490496131625075;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_22_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_22_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_22_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_22_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_22_var._parameters.mosaicity = 20;
  _Analyzer_22_var._parameters.domainthickness = 10;
  _Analyzer_22_var._parameters.temperature = 0;
  _Analyzer_22_var._parameters.verbose = 1;


  /* component Analyzer_22=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (39.71515800474511)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_22_var._rotation_absolute);
    rot_transpose(_Analyzer_21_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_22_var._rotation_absolute, tr1, _Analyzer_22_var._rotation_relative);
    _Analyzer_22_var._rotation_is_identity =  rot_test_identity(_Analyzer_22_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.800088458298231);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_22_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_21_var._position_absolute, _Analyzer_22_var._position_absolute);
    _Analyzer_22_var._position_relative = rot_apply(_Analyzer_22_var._rotation_absolute, tc1);
  } /* Analyzer_22=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_22", _Analyzer_22_var._position_absolute, _Analyzer_22_var._rotation_absolute);
  instrument->_position_absolute[27] = _Analyzer_22_var._position_absolute;
  instrument->_position_relative[27] = _Analyzer_22_var._position_relative;
    _Analyzer_22_var._position_relative_is_zero =  coords_test_zero(_Analyzer_22_var._position_relative);
  instrument->counter_N[27]  = instrument->counter_P[27] = instrument->counter_P2[27] = 0;
  instrument->counter_AbsorbProp[27]= 0;
  return(0);
} /* _Analyzer_22_setpos */

/* component Analyzer_23=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_23_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_23_setpos] component Analyzer_23=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_23_var._name, "Analyzer_23", 16384);
  stracpy(_Analyzer_23_var._type, "Monochromator_bent", 16384);
  _Analyzer_23_var._index=28;
  int current_setpos_index = 28;
  _Analyzer_23_var._parameters.zwidth = 0.08686616036731054;
  _Analyzer_23_var._parameters.yheight = 0.034140378549422014;
  _Analyzer_23_var._parameters.xthickness = 0.0005;
  _Analyzer_23_var._parameters.radius_x = 2.938196630867787;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_23_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_23_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_23_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_23_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_23_var._parameters.mosaicity = 20;
  _Analyzer_23_var._parameters.domainthickness = 10;
  _Analyzer_23_var._parameters.temperature = 0;
  _Analyzer_23_var._parameters.verbose = 1;


  /* component Analyzer_23=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (39.95707432499341)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_23_var._rotation_absolute);
    rot_transpose(_Analyzer_22_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_23_var._rotation_absolute, tr1, _Analyzer_23_var._rotation_relative);
    _Analyzer_23_var._rotation_is_identity =  rot_test_identity(_Analyzer_23_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8100884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_23_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_22_var._position_absolute, _Analyzer_23_var._position_absolute);
    _Analyzer_23_var._position_relative = rot_apply(_Analyzer_23_var._rotation_absolute, tc1);
  } /* Analyzer_23=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_23", _Analyzer_23_var._position_absolute, _Analyzer_23_var._rotation_absolute);
  instrument->_position_absolute[28] = _Analyzer_23_var._position_absolute;
  instrument->_position_relative[28] = _Analyzer_23_var._position_relative;
    _Analyzer_23_var._position_relative_is_zero =  coords_test_zero(_Analyzer_23_var._position_relative);
  instrument->counter_N[28]  = instrument->counter_P[28] = instrument->counter_P2[28] = 0;
  instrument->counter_AbsorbProp[28]= 0;
  return(0);
} /* _Analyzer_23_setpos */

/* component Analyzer_24=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_24_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_24_setpos] component Analyzer_24=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_24_var._name, "Analyzer_24", 16384);
  stracpy(_Analyzer_24_var._type, "Monochromator_bent", 16384);
  _Analyzer_24_var._index=29;
  int current_setpos_index = 29;
  _Analyzer_24_var._parameters.zwidth = 0.08644066416053885;
  _Analyzer_24_var._parameters.yheight = 0.03458712896971455;
  _Analyzer_24_var._parameters.xthickness = 0.0005;
  _Analyzer_24_var._parameters.radius_x = 2.9245083390640274;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_24_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_24_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_24_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_24_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_24_var._parameters.mosaicity = 20;
  _Analyzer_24_var._parameters.domainthickness = 10;
  _Analyzer_24_var._parameters.temperature = 0;
  _Analyzer_24_var._parameters.verbose = 1;


  /* component Analyzer_24=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (40.20037504749632)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_24_var._rotation_absolute);
    rot_transpose(_Analyzer_23_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_24_var._rotation_absolute, tr1, _Analyzer_24_var._rotation_relative);
    _Analyzer_24_var._rotation_is_identity =  rot_test_identity(_Analyzer_24_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8200884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_24_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_23_var._position_absolute, _Analyzer_24_var._position_absolute);
    _Analyzer_24_var._position_relative = rot_apply(_Analyzer_24_var._rotation_absolute, tc1);
  } /* Analyzer_24=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_24", _Analyzer_24_var._position_absolute, _Analyzer_24_var._rotation_absolute);
  instrument->_position_absolute[29] = _Analyzer_24_var._position_absolute;
  instrument->_position_relative[29] = _Analyzer_24_var._position_relative;
    _Analyzer_24_var._position_relative_is_zero =  coords_test_zero(_Analyzer_24_var._position_relative);
  instrument->counter_N[29]  = instrument->counter_P[29] = instrument->counter_P2[29] = 0;
  instrument->counter_AbsorbProp[29]= 0;
  return(0);
} /* _Analyzer_24_setpos */

/* component Analyzer_25=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_25_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_25_setpos] component Analyzer_25=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_25_var._name, "Analyzer_25", 16384);
  stracpy(_Analyzer_25_var._type, "Monochromator_bent", 16384);
  _Analyzer_25_var._index=30;
  int current_setpos_index = 30;
  _Analyzer_25_var._parameters.zwidth = 0.08742036977105623;
  _Analyzer_25_var._parameters.yheight = 0.03503379485569321;
  _Analyzer_25_var._parameters.xthickness = 0.0005;
  _Analyzer_25_var._parameters.radius_x = 2.898499909585128;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_25_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_25_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_25_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_25_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_25_var._parameters.mosaicity = 20;
  _Analyzer_25_var._parameters.domainthickness = 10;
  _Analyzer_25_var._parameters.temperature = 0;
  _Analyzer_25_var._parameters.verbose = 1;


  /* component Analyzer_25=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (40.44487724382992)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_25_var._rotation_absolute);
    rot_transpose(_Analyzer_24_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_25_var._rotation_absolute, tr1, _Analyzer_25_var._rotation_relative);
    _Analyzer_25_var._rotation_is_identity =  rot_test_identity(_Analyzer_25_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8300884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_25_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_24_var._position_absolute, _Analyzer_25_var._position_absolute);
    _Analyzer_25_var._position_relative = rot_apply(_Analyzer_25_var._rotation_absolute, tc1);
  } /* Analyzer_25=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_25", _Analyzer_25_var._position_absolute, _Analyzer_25_var._rotation_absolute);
  instrument->_position_absolute[30] = _Analyzer_25_var._position_absolute;
  instrument->_position_relative[30] = _Analyzer_25_var._position_relative;
    _Analyzer_25_var._position_relative_is_zero =  coords_test_zero(_Analyzer_25_var._position_relative);
  instrument->counter_N[30]  = instrument->counter_P[30] = instrument->counter_P2[30] = 0;
  instrument->counter_AbsorbProp[30]= 0;
  return(0);
} /* _Analyzer_25_setpos */

/* component Analyzer_26=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_26_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_26_setpos] component Analyzer_26=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_26_var._name, "Analyzer_26", 16384);
  stracpy(_Analyzer_26_var._type, "Monochromator_bent", 16384);
  _Analyzer_26_var._index=31;
  int current_setpos_index = 31;
  _Analyzer_26_var._parameters.zwidth = 0.08843868783512303;
  _Analyzer_26_var._parameters.yheight = 0.03543203625324905;
  _Analyzer_26_var._parameters.xthickness = 0.0005;
  _Analyzer_26_var._parameters.radius_x = 2.8712609724824634;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_26_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_26_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_26_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_26_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_26_var._parameters.mosaicity = 20;
  _Analyzer_26_var._parameters.domainthickness = 10;
  _Analyzer_26_var._parameters.temperature = 0;
  _Analyzer_26_var._parameters.verbose = 1;


  /* component Analyzer_26=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (40.69091649895137)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_26_var._rotation_absolute);
    rot_transpose(_Analyzer_25_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_26_var._rotation_absolute, tr1, _Analyzer_26_var._rotation_relative);
    _Analyzer_26_var._rotation_is_identity =  rot_test_identity(_Analyzer_26_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8400884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_26_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_25_var._position_absolute, _Analyzer_26_var._position_absolute);
    _Analyzer_26_var._position_relative = rot_apply(_Analyzer_26_var._rotation_absolute, tc1);
  } /* Analyzer_26=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_26", _Analyzer_26_var._position_absolute, _Analyzer_26_var._rotation_absolute);
  instrument->_position_absolute[31] = _Analyzer_26_var._position_absolute;
  instrument->_position_relative[31] = _Analyzer_26_var._position_relative;
    _Analyzer_26_var._position_relative_is_zero =  coords_test_zero(_Analyzer_26_var._position_relative);
  instrument->counter_N[31]  = instrument->counter_P[31] = instrument->counter_P2[31] = 0;
  instrument->counter_AbsorbProp[31]= 0;
  return(0);
} /* _Analyzer_26_setpos */

/* component Analyzer_27=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_27_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_27_setpos] component Analyzer_27=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_27_var._name, "Analyzer_27", 16384);
  stracpy(_Analyzer_27_var._type, "Monochromator_bent", 16384);
  _Analyzer_27_var._index=32;
  int current_setpos_index = 32;
  _Analyzer_27_var._parameters.zwidth = 0.08800627323583599;
  _Analyzer_27_var._parameters.yheight = 0.03587897867572419;
  _Analyzer_27_var._parameters.xthickness = 0.0005;
  _Analyzer_27_var._parameters.radius_x = 2.8357350040135745;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_27_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_27_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_27_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_27_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_27_var._parameters.mosaicity = 20;
  _Analyzer_27_var._parameters.domainthickness = 10;
  _Analyzer_27_var._parameters.temperature = 0;
  _Analyzer_27_var._parameters.verbose = 1;


  /* component Analyzer_27=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (40.938450294909146)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_27_var._rotation_absolute);
    rot_transpose(_Analyzer_26_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_27_var._rotation_absolute, tr1, _Analyzer_27_var._rotation_relative);
    _Analyzer_27_var._rotation_is_identity =  rot_test_identity(_Analyzer_27_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8500884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_27_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_26_var._position_absolute, _Analyzer_27_var._position_absolute);
    _Analyzer_27_var._position_relative = rot_apply(_Analyzer_27_var._rotation_absolute, tc1);
  } /* Analyzer_27=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_27", _Analyzer_27_var._position_absolute, _Analyzer_27_var._rotation_absolute);
  instrument->_position_absolute[32] = _Analyzer_27_var._position_absolute;
  instrument->_position_relative[32] = _Analyzer_27_var._position_relative;
    _Analyzer_27_var._position_relative_is_zero =  coords_test_zero(_Analyzer_27_var._position_relative);
  instrument->counter_N[32]  = instrument->counter_P[32] = instrument->counter_P2[32] = 0;
  instrument->counter_AbsorbProp[32]= 0;
  return(0);
} /* _Analyzer_27_setpos */

/* component Analyzer_28=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_28_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_28_setpos] component Analyzer_28=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_28_var._name, "Analyzer_28", 16384);
  stracpy(_Analyzer_28_var._type, "Monochromator_bent", 16384);
  _Analyzer_28_var._index=33;
  int current_setpos_index = 33;
  _Analyzer_28_var._parameters.zwidth = 0.08929978872070411;
  _Analyzer_28_var._parameters.yheight = 0.03632585465129224;
  _Analyzer_28_var._parameters.xthickness = 0.0005;
  _Analyzer_28_var._parameters.radius_x = 2.777987522899187;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_28_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_28_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_28_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_28_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_28_var._parameters.mosaicity = 20;
  _Analyzer_28_var._parameters.domainthickness = 10;
  _Analyzer_28_var._parameters.temperature = 0;
  _Analyzer_28_var._parameters.verbose = 1;


  /* component Analyzer_28=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (41.187885391194605)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_28_var._rotation_absolute);
    rot_transpose(_Analyzer_27_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_28_var._rotation_absolute, tr1, _Analyzer_28_var._rotation_relative);
    _Analyzer_28_var._rotation_is_identity =  rot_test_identity(_Analyzer_28_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8600884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_28_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_27_var._position_absolute, _Analyzer_28_var._position_absolute);
    _Analyzer_28_var._position_relative = rot_apply(_Analyzer_28_var._rotation_absolute, tc1);
  } /* Analyzer_28=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_28", _Analyzer_28_var._position_absolute, _Analyzer_28_var._rotation_absolute);
  instrument->_position_absolute[33] = _Analyzer_28_var._position_absolute;
  instrument->_position_relative[33] = _Analyzer_28_var._position_relative;
    _Analyzer_28_var._position_relative_is_zero =  coords_test_zero(_Analyzer_28_var._position_relative);
  instrument->counter_N[33]  = instrument->counter_P[33] = instrument->counter_P2[33] = 0;
  instrument->counter_AbsorbProp[33]= 0;
  return(0);
} /* _Analyzer_28_setpos */

/* component Analyzer_29=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_29_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_29_setpos] component Analyzer_29=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_29_var._name, "Analyzer_29", 16384);
  stracpy(_Analyzer_29_var._type, "Monochromator_bent", 16384);
  _Analyzer_29_var._index=34;
  int current_setpos_index = 34;
  _Analyzer_29_var._parameters.zwidth = 0.09028494113871582;
  _Analyzer_29_var._parameters.yheight = 0.036725534809827196;
  _Analyzer_29_var._parameters.xthickness = 0.0005;
  _Analyzer_29_var._parameters.radius_x = 2.728494211902119;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_29_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_29_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_29_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_29_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_29_var._parameters.mosaicity = 20;
  _Analyzer_29_var._parameters.domainthickness = 10;
  _Analyzer_29_var._parameters.temperature = 0;
  _Analyzer_29_var._parameters.verbose = 1;


  /* component Analyzer_29=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (41.44009035724288)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_29_var._rotation_absolute);
    rot_transpose(_Analyzer_28_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_29_var._rotation_absolute, tr1, _Analyzer_29_var._rotation_relative);
    _Analyzer_29_var._rotation_is_identity =  rot_test_identity(_Analyzer_29_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8700884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_29_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_28_var._position_absolute, _Analyzer_29_var._position_absolute);
    _Analyzer_29_var._position_relative = rot_apply(_Analyzer_29_var._rotation_absolute, tc1);
  } /* Analyzer_29=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_29", _Analyzer_29_var._position_absolute, _Analyzer_29_var._rotation_absolute);
  instrument->_position_absolute[34] = _Analyzer_29_var._position_absolute;
  instrument->_position_relative[34] = _Analyzer_29_var._position_relative;
    _Analyzer_29_var._position_relative_is_zero =  coords_test_zero(_Analyzer_29_var._position_relative);
  instrument->counter_N[34]  = instrument->counter_P[34] = instrument->counter_P2[34] = 0;
  instrument->counter_AbsorbProp[34]= 0;
  return(0);
} /* _Analyzer_29_setpos */

/* component Analyzer_30=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_30_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_30_setpos] component Analyzer_30=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_30_var._name, "Analyzer_30", 16384);
  stracpy(_Analyzer_30_var._type, "Monochromator_bent", 16384);
  _Analyzer_30_var._index=35;
  int current_setpos_index = 35;
  _Analyzer_30_var._parameters.zwidth = 0.08983821683994857;
  _Analyzer_30_var._parameters.yheight = 0.03717271464992166;
  _Analyzer_30_var._parameters.xthickness = 0.0005;
  _Analyzer_30_var._parameters.radius_x = 2.6637294000350398;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_30_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_30_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_30_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_30_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_30_var._parameters.mosaicity = 20;
  _Analyzer_30_var._parameters.domainthickness = 10;
  _Analyzer_30_var._parameters.temperature = 0;
  _Analyzer_30_var._parameters.verbose = 1;


  /* component Analyzer_30=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (41.69410145585204)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_30_var._rotation_absolute);
    rot_transpose(_Analyzer_29_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_30_var._rotation_absolute, tr1, _Analyzer_30_var._rotation_relative);
    _Analyzer_30_var._rotation_is_identity =  rot_test_identity(_Analyzer_30_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8800884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_30_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_29_var._position_absolute, _Analyzer_30_var._position_absolute);
    _Analyzer_30_var._position_relative = rot_apply(_Analyzer_30_var._rotation_absolute, tc1);
  } /* Analyzer_30=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_30", _Analyzer_30_var._position_absolute, _Analyzer_30_var._rotation_absolute);
  instrument->_position_absolute[35] = _Analyzer_30_var._position_absolute;
  instrument->_position_relative[35] = _Analyzer_30_var._position_relative;
    _Analyzer_30_var._position_relative_is_zero =  coords_test_zero(_Analyzer_30_var._position_relative);
  instrument->counter_N[35]  = instrument->counter_P[35] = instrument->counter_P2[35] = 0;
  instrument->counter_AbsorbProp[35]= 0;
  return(0);
} /* _Analyzer_30_setpos */

/* component Analyzer_31=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_31_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_31_setpos] component Analyzer_31=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_31_var._name, "Analyzer_31", 16384);
  stracpy(_Analyzer_31_var._type, "Monochromator_bent", 16384);
  _Analyzer_31_var._index=36;
  int current_setpos_index = 36;
  _Analyzer_31_var._parameters.zwidth = 0.09074428116517753;
  _Analyzer_31_var._parameters.yheight = 0.037619849175993314;
  _Analyzer_31_var._parameters.xthickness = 0.0005;
  _Analyzer_31_var._parameters.radius_x = 2.5682557576984952;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_31_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_31_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_31_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_31_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_31_var._parameters.mosaicity = 20;
  _Analyzer_31_var._parameters.domainthickness = 10;
  _Analyzer_31_var._parameters.temperature = 0;
  _Analyzer_31_var._parameters.verbose = 1;


  /* component Analyzer_31=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (41.951140797322154)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_31_var._rotation_absolute);
    rot_transpose(_Analyzer_30_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_31_var._rotation_absolute, tr1, _Analyzer_31_var._rotation_relative);
    _Analyzer_31_var._rotation_is_identity =  rot_test_identity(_Analyzer_31_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.8900884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_31_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_30_var._position_absolute, _Analyzer_31_var._position_absolute);
    _Analyzer_31_var._position_relative = rot_apply(_Analyzer_31_var._rotation_absolute, tc1);
  } /* Analyzer_31=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_31", _Analyzer_31_var._position_absolute, _Analyzer_31_var._rotation_absolute);
  instrument->_position_absolute[36] = _Analyzer_31_var._position_absolute;
  instrument->_position_relative[36] = _Analyzer_31_var._position_relative;
    _Analyzer_31_var._position_relative_is_zero =  coords_test_zero(_Analyzer_31_var._position_relative);
  instrument->counter_N[36]  = instrument->counter_P[36] = instrument->counter_P2[36] = 0;
  instrument->counter_AbsorbProp[36]= 0;
  return(0);
} /* _Analyzer_31_setpos */

/* component Analyzer_32=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_32_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_32_setpos] component Analyzer_32=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_32_var._name, "Analyzer_32", 16384);
  stracpy(_Analyzer_32_var._type, "Monochromator_bent", 16384);
  _Analyzer_32_var._index=37;
  int current_setpos_index = 37;
  _Analyzer_32_var._parameters.zwidth = 0.09203886640904534;
  _Analyzer_32_var._parameters.yheight = 0.03800950084407475;
  _Analyzer_32_var._parameters.xthickness = 0.0005;
  _Analyzer_32_var._parameters.radius_x = 2.4879223353464175;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_32_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_32_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_32_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_32_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_32_var._parameters.mosaicity = 20;
  _Analyzer_32_var._parameters.domainthickness = 10;
  _Analyzer_32_var._parameters.temperature = 0;
  _Analyzer_32_var._parameters.verbose = 1;


  /* component Analyzer_32=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (42.21348611029463)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_32_var._rotation_absolute);
    rot_transpose(_Analyzer_31_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_32_var._rotation_absolute, tr1, _Analyzer_32_var._rotation_relative);
    _Analyzer_32_var._rotation_is_identity =  rot_test_identity(_Analyzer_32_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9000884582982311);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_32_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_31_var._position_absolute, _Analyzer_32_var._position_absolute);
    _Analyzer_32_var._position_relative = rot_apply(_Analyzer_32_var._rotation_absolute, tc1);
  } /* Analyzer_32=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_32", _Analyzer_32_var._position_absolute, _Analyzer_32_var._rotation_absolute);
  instrument->_position_absolute[37] = _Analyzer_32_var._position_absolute;
  instrument->_position_relative[37] = _Analyzer_32_var._position_relative;
    _Analyzer_32_var._position_relative_is_zero =  coords_test_zero(_Analyzer_32_var._position_relative);
  instrument->counter_N[37]  = instrument->counter_P[37] = instrument->counter_P2[37] = 0;
  instrument->counter_AbsorbProp[37]= 0;
  return(0);
} /* _Analyzer_32_setpos */

/* component Analyzer_33=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_33_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_33_setpos] component Analyzer_33=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_33_var._name, "Analyzer_33", 16384);
  stracpy(_Analyzer_33_var._type, "Monochromator_bent", 16384);
  _Analyzer_33_var._index=38;
  int current_setpos_index = 38;
  _Analyzer_33_var._parameters.zwidth = 0.092558807854864;
  _Analyzer_33_var._parameters.yheight = 0.03845708019977417;
  _Analyzer_33_var._parameters.xthickness = 0.0005;
  _Analyzer_33_var._parameters.radius_x = 2.3557164020632526;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_33_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_33_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_33_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_33_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_33_var._parameters.mosaicity = 20;
  _Analyzer_33_var._parameters.domainthickness = 10;
  _Analyzer_33_var._parameters.temperature = 0;
  _Analyzer_33_var._parameters.verbose = 1;


  /* component Analyzer_33=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (42.47873655390719)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_33_var._rotation_absolute);
    rot_transpose(_Analyzer_32_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_33_var._rotation_absolute, tr1, _Analyzer_33_var._rotation_relative);
    _Analyzer_33_var._rotation_is_identity =  rot_test_identity(_Analyzer_33_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9100884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_33_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_32_var._position_absolute, _Analyzer_33_var._position_absolute);
    _Analyzer_33_var._position_relative = rot_apply(_Analyzer_33_var._rotation_absolute, tc1);
  } /* Analyzer_33=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_33", _Analyzer_33_var._position_absolute, _Analyzer_33_var._rotation_absolute);
  instrument->_position_absolute[38] = _Analyzer_33_var._position_absolute;
  instrument->_position_relative[38] = _Analyzer_33_var._position_relative;
    _Analyzer_33_var._position_relative_is_zero =  coords_test_zero(_Analyzer_33_var._position_relative);
  instrument->counter_N[38]  = instrument->counter_P[38] = instrument->counter_P2[38] = 0;
  instrument->counter_AbsorbProp[38]= 0;
  return(0);
} /* _Analyzer_33_setpos */

/* component Analyzer_34=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_34_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_34_setpos] component Analyzer_34=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_34_var._name, "Analyzer_34", 16384);
  stracpy(_Analyzer_34_var._type, "Monochromator_bent", 16384);
  _Analyzer_34_var._index=39;
  int current_setpos_index = 39;
  _Analyzer_34_var._parameters.zwidth = 0.08809032316180643;
  _Analyzer_34_var._parameters.yheight = 0.03890464199783919;
  _Analyzer_34_var._parameters.xthickness = 0.0005;
  _Analyzer_34_var._parameters.radius_x = 2.3145921881417686;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_34_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_34_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_34_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_34_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_34_var._parameters.mosaicity = 20;
  _Analyzer_34_var._parameters.domainthickness = 10;
  _Analyzer_34_var._parameters.temperature = 0;
  _Analyzer_34_var._parameters.verbose = 1;


  /* component Analyzer_34=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (42.746219505063564)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_34_var._rotation_absolute);
    rot_transpose(_Analyzer_33_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_34_var._rotation_absolute, tr1, _Analyzer_34_var._rotation_relative);
    _Analyzer_34_var._rotation_is_identity =  rot_test_identity(_Analyzer_34_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9200884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_34_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_33_var._position_absolute, _Analyzer_34_var._position_absolute);
    _Analyzer_34_var._position_relative = rot_apply(_Analyzer_34_var._rotation_absolute, tc1);
  } /* Analyzer_34=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_34", _Analyzer_34_var._position_absolute, _Analyzer_34_var._rotation_absolute);
  instrument->_position_absolute[39] = _Analyzer_34_var._position_absolute;
  instrument->_position_relative[39] = _Analyzer_34_var._position_relative;
    _Analyzer_34_var._position_relative_is_zero =  coords_test_zero(_Analyzer_34_var._position_relative);
  instrument->counter_N[39]  = instrument->counter_P[39] = instrument->counter_P2[39] = 0;
  instrument->counter_AbsorbProp[39]= 0;
  return(0);
} /* _Analyzer_34_setpos */

/* component Analyzer_35=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_35_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_35_setpos] component Analyzer_35=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_35_var._name, "Analyzer_35", 16384);
  stracpy(_Analyzer_35_var._type, "Monochromator_bent", 16384);
  _Analyzer_35_var._index=40;
  int current_setpos_index = 40;
  _Analyzer_35_var._parameters.zwidth = 0.08470875813671634;
  _Analyzer_35_var._parameters.yheight = 0.03930740330667787;
  _Analyzer_35_var._parameters.xthickness = 0.0005;
  _Analyzer_35_var._parameters.radius_x = 2.291483856097038;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_35_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_35_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_35_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_35_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_35_var._parameters.mosaicity = 20;
  _Analyzer_35_var._parameters.domainthickness = 10;
  _Analyzer_35_var._parameters.temperature = 0;
  _Analyzer_35_var._parameters.verbose = 1;


  /* component Analyzer_35=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (43.01875584443871)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_35_var._rotation_absolute);
    rot_transpose(_Analyzer_34_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_35_var._rotation_absolute, tr1, _Analyzer_35_var._rotation_relative);
    _Analyzer_35_var._rotation_is_identity =  rot_test_identity(_Analyzer_35_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9300884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_35_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_34_var._position_absolute, _Analyzer_35_var._position_absolute);
    _Analyzer_35_var._position_relative = rot_apply(_Analyzer_35_var._rotation_absolute, tc1);
  } /* Analyzer_35=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_35", _Analyzer_35_var._position_absolute, _Analyzer_35_var._rotation_absolute);
  instrument->_position_absolute[40] = _Analyzer_35_var._position_absolute;
  instrument->_position_relative[40] = _Analyzer_35_var._position_relative;
    _Analyzer_35_var._position_relative_is_zero =  coords_test_zero(_Analyzer_35_var._position_relative);
  instrument->counter_N[40]  = instrument->counter_P[40] = instrument->counter_P2[40] = 0;
  instrument->counter_AbsorbProp[40]= 0;
  return(0);
} /* _Analyzer_35_setpos */

/* component Analyzer_36=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_36_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_36_setpos] component Analyzer_36=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_36_var._name, "Analyzer_36", 16384);
  stracpy(_Analyzer_36_var._type, "Monochromator_bent", 16384);
  _Analyzer_36_var._index=41;
  int current_setpos_index = 41;
  _Analyzer_36_var._parameters.zwidth = 0.0799700742586517;
  _Analyzer_36_var._parameters.yheight = 0.03975537102141946;
  _Analyzer_36_var._parameters.xthickness = 0.0005;
  _Analyzer_36_var._parameters.radius_x = 2.2505408566585348;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_36_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_36_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_36_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_36_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_36_var._parameters.mosaicity = 20;
  _Analyzer_36_var._parameters.domainthickness = 10;
  _Analyzer_36_var._parameters.temperature = 0;
  _Analyzer_36_var._parameters.verbose = 1;


  /* component Analyzer_36=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (43.2930526461892)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_36_var._rotation_absolute);
    rot_transpose(_Analyzer_35_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_36_var._rotation_absolute, tr1, _Analyzer_36_var._rotation_relative);
    _Analyzer_36_var._rotation_is_identity =  rot_test_identity(_Analyzer_36_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9400884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_36_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_35_var._position_absolute, _Analyzer_36_var._position_absolute);
    _Analyzer_36_var._position_relative = rot_apply(_Analyzer_36_var._rotation_absolute, tc1);
  } /* Analyzer_36=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_36", _Analyzer_36_var._position_absolute, _Analyzer_36_var._rotation_absolute);
  instrument->_position_absolute[41] = _Analyzer_36_var._position_absolute;
  instrument->_position_relative[41] = _Analyzer_36_var._position_relative;
    _Analyzer_36_var._position_relative_is_zero =  coords_test_zero(_Analyzer_36_var._position_relative);
  instrument->counter_N[41]  = instrument->counter_P[41] = instrument->counter_P2[41] = 0;
  instrument->counter_AbsorbProp[41]= 0;
  return(0);
} /* _Analyzer_36_setpos */

/* component Analyzer_37=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_37_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_37_setpos] component Analyzer_37=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_37_var._name, "Analyzer_37", 16384);
  stracpy(_Analyzer_37_var._type, "Monochromator_bent", 16384);
  _Analyzer_37_var._index=42;
  int current_setpos_index = 42;
  _Analyzer_37_var._parameters.zwidth = 0.07793790704996095;
  _Analyzer_37_var._parameters.yheight = 0.040170340093780246;
  _Analyzer_37_var._parameters.xthickness = 0.0005;
  _Analyzer_37_var._parameters.radius_x = 2.1756847577357843;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_37_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_37_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_37_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_37_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_37_var._parameters.mosaicity = 20;
  _Analyzer_37_var._parameters.domainthickness = 10;
  _Analyzer_37_var._parameters.temperature = 0;
  _Analyzer_37_var._parameters.verbose = 1;


  /* component Analyzer_37=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (43.577221074552455)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_37_var._rotation_absolute);
    rot_transpose(_Analyzer_36_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_37_var._rotation_absolute, tr1, _Analyzer_37_var._rotation_relative);
    _Analyzer_37_var._rotation_is_identity =  rot_test_identity(_Analyzer_37_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9500884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_37_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_36_var._position_absolute, _Analyzer_37_var._position_absolute);
    _Analyzer_37_var._position_relative = rot_apply(_Analyzer_37_var._rotation_absolute, tc1);
  } /* Analyzer_37=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_37", _Analyzer_37_var._position_absolute, _Analyzer_37_var._rotation_absolute);
  instrument->_position_absolute[42] = _Analyzer_37_var._position_absolute;
  instrument->_position_relative[42] = _Analyzer_37_var._position_relative;
    _Analyzer_37_var._position_relative_is_zero =  coords_test_zero(_Analyzer_37_var._position_relative);
  instrument->counter_N[42]  = instrument->counter_P[42] = instrument->counter_P2[42] = 0;
  instrument->counter_AbsorbProp[42]= 0;
  return(0);
} /* _Analyzer_37_setpos */

/* component Analyzer_38=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_38_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_38_setpos] component Analyzer_38=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_38_var._name, "Analyzer_38", 16384);
  stracpy(_Analyzer_38_var._type, "Monochromator_bent", 16384);
  _Analyzer_38_var._index=43;
  int current_setpos_index = 43;
  _Analyzer_38_var._parameters.zwidth = 0.07360030790545345;
  _Analyzer_38_var._parameters.yheight = 0.04060776698549516;
  _Analyzer_38_var._parameters.xthickness = 0.0005;
  _Analyzer_38_var._parameters.radius_x = 2.1295483923353666;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_38_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_38_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_38_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_38_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_38_var._parameters.mosaicity = 20;
  _Analyzer_38_var._parameters.domainthickness = 10;
  _Analyzer_38_var._parameters.temperature = 0;
  _Analyzer_38_var._parameters.verbose = 1;


  /* component Analyzer_38=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (43.8620488629705)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_38_var._rotation_absolute);
    rot_transpose(_Analyzer_37_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_38_var._rotation_absolute, tr1, _Analyzer_38_var._rotation_relative);
    _Analyzer_38_var._rotation_is_identity =  rot_test_identity(_Analyzer_38_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9600884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_38_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_37_var._position_absolute, _Analyzer_38_var._position_absolute);
    _Analyzer_38_var._position_relative = rot_apply(_Analyzer_38_var._rotation_absolute, tc1);
  } /* Analyzer_38=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_38", _Analyzer_38_var._position_absolute, _Analyzer_38_var._rotation_absolute);
  instrument->_position_absolute[43] = _Analyzer_38_var._position_absolute;
  instrument->_position_relative[43] = _Analyzer_38_var._position_relative;
    _Analyzer_38_var._position_relative_is_zero =  coords_test_zero(_Analyzer_38_var._position_relative);
  instrument->counter_N[43]  = instrument->counter_P[43] = instrument->counter_P2[43] = 0;
  instrument->counter_AbsorbProp[43]= 0;
  return(0);
} /* _Analyzer_38_setpos */

/* component Analyzer_39=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_39_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_39_setpos] component Analyzer_39=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_39_var._name, "Analyzer_39", 16384);
  stracpy(_Analyzer_39_var._type, "Monochromator_bent", 16384);
  _Analyzer_39_var._index=44;
  int current_setpos_index = 44;
  _Analyzer_39_var._parameters.zwidth = 0.06895356082425434;
  _Analyzer_39_var._parameters.yheight = 0.04105628146501619;
  _Analyzer_39_var._parameters.xthickness = 0.0005;
  _Analyzer_39_var._parameters.radius_x = 2.0741260312112377;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_39_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_39_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_39_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_39_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_39_var._parameters.mosaicity = 20;
  _Analyzer_39_var._parameters.domainthickness = 10;
  _Analyzer_39_var._parameters.temperature = 0;
  _Analyzer_39_var._parameters.verbose = 1;


  /* component Analyzer_39=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (44.15256729065363)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_39_var._rotation_absolute);
    rot_transpose(_Analyzer_38_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_39_var._rotation_absolute, tr1, _Analyzer_39_var._rotation_relative);
    _Analyzer_39_var._rotation_is_identity =  rot_test_identity(_Analyzer_39_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9700884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_39_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_38_var._position_absolute, _Analyzer_39_var._position_absolute);
    _Analyzer_39_var._position_relative = rot_apply(_Analyzer_39_var._rotation_absolute, tc1);
  } /* Analyzer_39=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_39", _Analyzer_39_var._position_absolute, _Analyzer_39_var._rotation_absolute);
  instrument->_position_absolute[44] = _Analyzer_39_var._position_absolute;
  instrument->_position_relative[44] = _Analyzer_39_var._position_relative;
    _Analyzer_39_var._position_relative_is_zero =  coords_test_zero(_Analyzer_39_var._position_relative);
  instrument->counter_N[44]  = instrument->counter_P[44] = instrument->counter_P2[44] = 0;
  instrument->counter_AbsorbProp[44]= 0;
  return(0);
} /* _Analyzer_39_setpos */

/* component Analyzer_40=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_40_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_40_setpos] component Analyzer_40=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_40_var._name, "Analyzer_40", 16384);
  stracpy(_Analyzer_40_var._type, "Monochromator_bent", 16384);
  _Analyzer_40_var._index=45;
  int current_setpos_index = 45;
  _Analyzer_40_var._parameters.zwidth = 0.06603219877163607;
  _Analyzer_40_var._parameters.yheight = 0.04145128904284828;
  _Analyzer_40_var._parameters.xthickness = 0.0005;
  _Analyzer_40_var._parameters.radius_x = 2.036498894266056;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_40_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_40_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_40_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_40_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_40_var._parameters.mosaicity = 20;
  _Analyzer_40_var._parameters.domainthickness = 10;
  _Analyzer_40_var._parameters.temperature = 0;
  _Analyzer_40_var._parameters.verbose = 1;


  /* component Analyzer_40=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (44.449534626673206)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_40_var._rotation_absolute);
    rot_transpose(_Analyzer_39_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_40_var._rotation_absolute, tr1, _Analyzer_40_var._rotation_relative);
    _Analyzer_40_var._rotation_is_identity =  rot_test_identity(_Analyzer_40_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9800884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_40_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_39_var._position_absolute, _Analyzer_40_var._position_absolute);
    _Analyzer_40_var._position_relative = rot_apply(_Analyzer_40_var._rotation_absolute, tc1);
  } /* Analyzer_40=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_40", _Analyzer_40_var._position_absolute, _Analyzer_40_var._rotation_absolute);
  instrument->_position_absolute[45] = _Analyzer_40_var._position_absolute;
  instrument->_position_relative[45] = _Analyzer_40_var._position_relative;
    _Analyzer_40_var._position_relative_is_zero =  coords_test_zero(_Analyzer_40_var._position_relative);
  instrument->counter_N[45]  = instrument->counter_P[45] = instrument->counter_P2[45] = 0;
  instrument->counter_AbsorbProp[45]= 0;
  return(0);
} /* _Analyzer_40_setpos */

/* component Analyzer_41=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_41_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_41_setpos] component Analyzer_41=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_41_var._name, "Analyzer_41", 16384);
  stracpy(_Analyzer_41_var._type, "Monochromator_bent", 16384);
  _Analyzer_41_var._index=46;
  int current_setpos_index = 46;
  _Analyzer_41_var._parameters.zwidth = 0.061445022067135514;
  _Analyzer_41_var._parameters.yheight = 0.041900561961020134;
  _Analyzer_41_var._parameters.xthickness = 0.0005;
  _Analyzer_41_var._parameters.radius_x = 1.9705411161560504;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_41_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_41_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_41_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_41_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_41_var._parameters.mosaicity = 20;
  _Analyzer_41_var._parameters.domainthickness = 10;
  _Analyzer_41_var._parameters.temperature = 0;
  _Analyzer_41_var._parameters.verbose = 1;


  /* component Analyzer_41=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (44.754531942864546)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_41_var._rotation_absolute);
    rot_transpose(_Analyzer_40_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_41_var._rotation_absolute, tr1, _Analyzer_41_var._rotation_relative);
    _Analyzer_41_var._rotation_is_identity =  rot_test_identity(_Analyzer_41_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9900884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_41_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_40_var._position_absolute, _Analyzer_41_var._position_absolute);
    _Analyzer_41_var._position_relative = rot_apply(_Analyzer_41_var._rotation_absolute, tc1);
  } /* Analyzer_41=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_41", _Analyzer_41_var._position_absolute, _Analyzer_41_var._rotation_absolute);
  instrument->_position_absolute[46] = _Analyzer_41_var._position_absolute;
  instrument->_position_relative[46] = _Analyzer_41_var._position_relative;
    _Analyzer_41_var._position_relative_is_zero =  coords_test_zero(_Analyzer_41_var._position_relative);
  instrument->counter_N[46]  = instrument->counter_P[46] = instrument->counter_P2[46] = 0;
  instrument->counter_AbsorbProp[46]= 0;
  return(0);
} /* _Analyzer_41_setpos */

/* component Analyzer_42=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_42_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_42_setpos] component Analyzer_42=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_42_var._name, "Analyzer_42", 16384);
  stracpy(_Analyzer_42_var._type, "Monochromator_bent", 16384);
  _Analyzer_42_var._index=47;
  int current_setpos_index = 47;
  _Analyzer_42_var._parameters.zwidth = 0.056891594444929736;
  _Analyzer_42_var._parameters.yheight = 0.04235005124793069;
  _Analyzer_42_var._parameters.xthickness = 0.0005;
  _Analyzer_42_var._parameters.radius_x = 1.897817106671356;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_42_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_42_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_42_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_42_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_42_var._parameters.mosaicity = 20;
  _Analyzer_42_var._parameters.domainthickness = 10;
  _Analyzer_42_var._parameters.temperature = 0;
  _Analyzer_42_var._parameters.verbose = 1;


  /* component Analyzer_42=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (45.06947238563547)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_42_var._rotation_absolute);
    rot_transpose(_Analyzer_41_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_42_var._rotation_absolute, tr1, _Analyzer_42_var._rotation_relative);
    _Analyzer_42_var._rotation_is_identity =  rot_test_identity(_Analyzer_42_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1.0000884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_42_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_41_var._position_absolute, _Analyzer_42_var._position_absolute);
    _Analyzer_42_var._position_relative = rot_apply(_Analyzer_42_var._rotation_absolute, tc1);
  } /* Analyzer_42=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_42", _Analyzer_42_var._position_absolute, _Analyzer_42_var._rotation_absolute);
  instrument->_position_absolute[47] = _Analyzer_42_var._position_absolute;
  instrument->_position_relative[47] = _Analyzer_42_var._position_relative;
    _Analyzer_42_var._position_relative_is_zero =  coords_test_zero(_Analyzer_42_var._position_relative);
  instrument->counter_N[47]  = instrument->counter_P[47] = instrument->counter_P2[47] = 0;
  instrument->counter_AbsorbProp[47]= 0;
  return(0);
} /* _Analyzer_42_setpos */

/* component Analyzer_43=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_43_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_43_setpos] component Analyzer_43=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_43_var._name, "Analyzer_43", 16384);
  stracpy(_Analyzer_43_var._type, "Monochromator_bent", 16384);
  _Analyzer_43_var._index=48;
  int current_setpos_index = 48;
  _Analyzer_43_var._parameters.zwidth = 0.0537062314545185;
  _Analyzer_43_var._parameters.yheight = 0.04275810889984801;
  _Analyzer_43_var._parameters.xthickness = 0.0005;
  _Analyzer_43_var._parameters.radius_x = 1.8394584169039285;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_43_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_43_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_43_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_43_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_43_var._parameters.mosaicity = 20;
  _Analyzer_43_var._parameters.domainthickness = 10;
  _Analyzer_43_var._parameters.temperature = 0;
  _Analyzer_43_var._parameters.verbose = 1;


  /* component Analyzer_43=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (45.39286461363646)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_43_var._rotation_absolute);
    rot_transpose(_Analyzer_42_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_43_var._rotation_absolute, tr1, _Analyzer_43_var._rotation_relative);
    _Analyzer_43_var._rotation_is_identity =  rot_test_identity(_Analyzer_43_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1.0100884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_43_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_42_var._position_absolute, _Analyzer_43_var._position_absolute);
    _Analyzer_43_var._position_relative = rot_apply(_Analyzer_43_var._rotation_absolute, tc1);
  } /* Analyzer_43=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_43", _Analyzer_43_var._position_absolute, _Analyzer_43_var._rotation_absolute);
  instrument->_position_absolute[48] = _Analyzer_43_var._position_absolute;
  instrument->_position_relative[48] = _Analyzer_43_var._position_relative;
    _Analyzer_43_var._position_relative_is_zero =  coords_test_zero(_Analyzer_43_var._position_relative);
  instrument->counter_N[48]  = instrument->counter_P[48] = instrument->counter_P2[48] = 0;
  instrument->counter_AbsorbProp[48]= 0;
  return(0);
} /* _Analyzer_43_setpos */

/* component Analyzer_44=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_44_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_44_setpos] component Analyzer_44=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_44_var._name, "Analyzer_44", 16384);
  stracpy(_Analyzer_44_var._type, "Monochromator_bent", 16384);
  _Analyzer_44_var._index=49;
  int current_setpos_index = 49;
  _Analyzer_44_var._parameters.zwidth = 0.04824834288562266;
  _Analyzer_44_var._parameters.yheight = 0.04320877236345719;
  _Analyzer_44_var._parameters.xthickness = 0.0005;
  _Analyzer_44_var._parameters.radius_x = 1.784778344731803;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_44_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_44_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_44_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_44_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_44_var._parameters.mosaicity = 20;
  _Analyzer_44_var._parameters.domainthickness = 10;
  _Analyzer_44_var._parameters.temperature = 0;
  _Analyzer_44_var._parameters.verbose = 1;


  /* component Analyzer_44=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (45.719934361579654)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_44_var._rotation_absolute);
    rot_transpose(_Analyzer_43_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_44_var._rotation_absolute, tr1, _Analyzer_44_var._rotation_relative);
    _Analyzer_44_var._rotation_is_identity =  rot_test_identity(_Analyzer_44_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1.0200884582982312);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_44_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_43_var._position_absolute, _Analyzer_44_var._position_absolute);
    _Analyzer_44_var._position_relative = rot_apply(_Analyzer_44_var._rotation_absolute, tc1);
  } /* Analyzer_44=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_44", _Analyzer_44_var._position_absolute, _Analyzer_44_var._rotation_absolute);
  instrument->_position_absolute[49] = _Analyzer_44_var._position_absolute;
  instrument->_position_relative[49] = _Analyzer_44_var._position_relative;
    _Analyzer_44_var._position_relative_is_zero =  coords_test_zero(_Analyzer_44_var._position_relative);
  instrument->counter_N[49]  = instrument->counter_P[49] = instrument->counter_P2[49] = 0;
  instrument->counter_AbsorbProp[49]= 0;
  return(0);
} /* _Analyzer_44_setpos */

/* component Analyzer_45=Monochromator_bent() SETTING, POSITION/ROTATION */
int _Analyzer_45_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Analyzer_45_setpos] component Analyzer_45=Monochromator_bent() SETTING [Monochromator_bent.comp:739]");
  stracpy(_Analyzer_45_var._name, "Analyzer_45", 16384);
  stracpy(_Analyzer_45_var._type, "Monochromator_bent", 16384);
  _Analyzer_45_var._index=50;
  int current_setpos_index = 50;
  _Analyzer_45_var._parameters.zwidth = 0.04347135501377188;
  _Analyzer_45_var._parameters.yheight = 0.04366020910051233;
  _Analyzer_45_var._parameters.xthickness = 0.0005;
  _Analyzer_45_var._parameters.radius_x = 1.7009660804373992;
  if("Si111" && strlen("Si111"))
    stracpy(_Analyzer_45_var._parameters.plane_of_reflection, "Si111" ? "Si111" : "", 16384);
  else 
  _Analyzer_45_var._parameters.plane_of_reflection[0]='\0';
  _Analyzer_45_var._parameters.angle_to_cut_horizontal = 0;
  _Analyzer_45_var._parameters.angle_to_cut_vertical = 0;
  _Analyzer_45_var._parameters.mosaicity = 20;
  _Analyzer_45_var._parameters.domainthickness = 10;
  _Analyzer_45_var._parameters.temperature = 0;
  _Analyzer_45_var._parameters.verbose = 1;


  /* component Analyzer_45=Monochromator_bent() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (46.070847274540796)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _Analyzer_45_var._rotation_absolute);
    rot_transpose(_Analyzer_44_var._rotation_absolute, tr1);
    rot_mul(_Analyzer_45_var._rotation_absolute, tr1, _Analyzer_45_var._rotation_relative);
    _Analyzer_45_var._rotation_is_identity =  rot_test_identity(_Analyzer_45_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1.0300884582982313);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Analyzer_45_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_44_var._position_absolute, _Analyzer_45_var._position_absolute);
    _Analyzer_45_var._position_relative = rot_apply(_Analyzer_45_var._rotation_absolute, tc1);
  } /* Analyzer_45=Monochromator_bent() AT ROTATED */
  DEBUG_COMPONENT("Analyzer_45", _Analyzer_45_var._position_absolute, _Analyzer_45_var._rotation_absolute);
  instrument->_position_absolute[50] = _Analyzer_45_var._position_absolute;
  instrument->_position_relative[50] = _Analyzer_45_var._position_relative;
    _Analyzer_45_var._position_relative_is_zero =  coords_test_zero(_Analyzer_45_var._position_relative);
  instrument->counter_N[50]  = instrument->counter_P[50] = instrument->counter_P2[50] = 0;
  instrument->counter_AbsorbProp[50]= 0;
  return(0);
} /* _Analyzer_45_setpos */

/* component rotator=Arm() SETTING, POSITION/ROTATION */
int _rotator_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_rotator_setpos] component rotator=Arm() SETTING [Arm:0]");
  stracpy(_rotator_var._name, "rotator", 16384);
  stracpy(_rotator_var._type, "Arm", 16384);
  _rotator_var._index=51;
  int current_setpos_index = 51;
  /* component rotator=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (-90)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _rotator_var._rotation_absolute);
    rot_transpose(_Analyzer_45_var._rotation_absolute, tr1);
    rot_mul(_rotator_var._rotation_absolute, tr1, _rotator_var._rotation_relative);
    _rotator_var._rotation_is_identity =  rot_test_identity(_rotator_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _rotator_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_45_var._position_absolute, _rotator_var._position_absolute);
    _rotator_var._position_relative = rot_apply(_rotator_var._rotation_absolute, tc1);
  } /* rotator=Arm() AT ROTATED */
  DEBUG_COMPONENT("rotator", _rotator_var._position_absolute, _rotator_var._rotation_absolute);
  instrument->_position_absolute[51] = _rotator_var._position_absolute;
  instrument->_position_relative[51] = _rotator_var._position_relative;
    _rotator_var._position_relative_is_zero =  coords_test_zero(_rotator_var._position_relative);
  instrument->counter_N[51]  = instrument->counter_P[51] = instrument->counter_P2[51] = 0;
  instrument->counter_AbsorbProp[51]= 0;
  return(0);
} /* _rotator_setpos */

/* component detector_pos=Arm() SETTING, POSITION/ROTATION */
int _detector_pos_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_detector_pos_setpos] component detector_pos=Arm() SETTING [Arm:0]");
  stracpy(_detector_pos_var._name, "detector_pos", 16384);
  stracpy(_detector_pos_var._type, "Arm", 16384);
  _detector_pos_var._index=52;
  int current_setpos_index = 52;
  /* component detector_pos=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _rotator_var._rotation_absolute, _detector_pos_var._rotation_absolute);
    rot_transpose(_Analyzer_45_var._rotation_absolute, tr1);
    rot_mul(_detector_pos_var._rotation_absolute, tr1, _detector_pos_var._rotation_relative);
    _detector_pos_var._rotation_is_identity =  rot_test_identity(_detector_pos_var._rotation_relative);
    tc1 = coords_set(
      0, _instrument_var._parameters.L1, _instrument_var._parameters.L0);
    rot_transpose(_rotator_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _detector_pos_var._position_absolute = coords_add(_rotator_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_45_var._position_absolute, _detector_pos_var._position_absolute);
    _detector_pos_var._position_relative = rot_apply(_detector_pos_var._rotation_absolute, tc1);
  } /* detector_pos=Arm() AT ROTATED */
  DEBUG_COMPONENT("detector_pos", _detector_pos_var._position_absolute, _detector_pos_var._rotation_absolute);
  instrument->_position_absolute[52] = _detector_pos_var._position_absolute;
  instrument->_position_relative[52] = _detector_pos_var._position_relative;
    _detector_pos_var._position_relative_is_zero =  coords_test_zero(_detector_pos_var._position_relative);
  instrument->counter_N[52]  = instrument->counter_P[52] = instrument->counter_P2[52] = 0;
  instrument->counter_AbsorbProp[52]= 0;
  return(0);
} /* _detector_pos_setpos */

/* component detector=Arm() SETTING, POSITION/ROTATION */
int _detector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_detector_setpos] component detector=Arm() SETTING [Arm:0]");
  stracpy(_detector_var._name, "detector", 16384);
  stracpy(_detector_var._type, "Arm", 16384);
  _detector_var._index=53;
  int current_setpos_index = 53;
  /* component detector=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (90 + _instrument_var._parameters.det_rot)*DEG2RAD, (0)*DEG2RAD, (90)*DEG2RAD);
    rot_mul(tr1, _detector_pos_var._rotation_absolute, _detector_var._rotation_absolute);
    rot_transpose(_Analyzer_45_var._rotation_absolute, tr1);
    rot_mul(_detector_var._rotation_absolute, tr1, _detector_var._rotation_relative);
    _detector_var._rotation_is_identity =  rot_test_identity(_detector_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_detector_pos_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _detector_var._position_absolute = coords_add(_detector_pos_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_45_var._position_absolute, _detector_var._position_absolute);
    _detector_var._position_relative = rot_apply(_detector_var._rotation_absolute, tc1);
  } /* detector=Arm() AT ROTATED */
  DEBUG_COMPONENT("detector", _detector_var._position_absolute, _detector_var._rotation_absolute);
  instrument->_position_absolute[53] = _detector_var._position_absolute;
  instrument->_position_relative[53] = _detector_var._position_relative;
    _detector_var._position_relative_is_zero =  coords_test_zero(_detector_var._position_relative);
  instrument->counter_N[53]  = instrument->counter_P[53] = instrument->counter_P2[53] = 0;
  instrument->counter_AbsorbProp[53]= 0;
  return(0);
} /* _detector_setpos */

/* component psd_monitor_end=PSD_monitor() SETTING, POSITION/ROTATION */
int _psd_monitor_end_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psd_monitor_end_setpos] component psd_monitor_end=PSD_monitor() SETTING [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:62]");
  stracpy(_psd_monitor_end_var._name, "psd_monitor_end", 16384);
  stracpy(_psd_monitor_end_var._type, "PSD_monitor", 16384);
  _psd_monitor_end_var._index=54;
  int current_setpos_index = 54;
  _psd_monitor_end_var._parameters.nx = 1001;
  _psd_monitor_end_var._parameters.ny = 1001;
  _psd_monitor_end_var._parameters.filename[0]='\0';
  _psd_monitor_end_var._parameters.xmin = -0.05;
  _psd_monitor_end_var._parameters.xmax = 0.05;
  _psd_monitor_end_var._parameters.ymin = -0.05;
  _psd_monitor_end_var._parameters.ymax = 0.05;
  _psd_monitor_end_var._parameters.xwidth = _instrument_var._parameters.Ld;
  _psd_monitor_end_var._parameters.yheight = 0.5;
  _psd_monitor_end_var._parameters.restore_neutron = 1;
  _psd_monitor_end_var._parameters.nowritefile = 0;


  /* component psd_monitor_end=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _detector_var._rotation_absolute, _psd_monitor_end_var._rotation_absolute);
    rot_transpose(_Analyzer_45_var._rotation_absolute, tr1);
    rot_mul(_psd_monitor_end_var._rotation_absolute, tr1, _psd_monitor_end_var._rotation_relative);
    _psd_monitor_end_var._rotation_is_identity =  rot_test_identity(_psd_monitor_end_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_detector_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psd_monitor_end_var._position_absolute = coords_add(_detector_var._position_absolute, tc2);
    tc1 = coords_sub(_Analyzer_45_var._position_absolute, _psd_monitor_end_var._position_absolute);
    _psd_monitor_end_var._position_relative = rot_apply(_psd_monitor_end_var._rotation_absolute, tc1);
  } /* psd_monitor_end=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psd_monitor_end", _psd_monitor_end_var._position_absolute, _psd_monitor_end_var._rotation_absolute);
  instrument->_position_absolute[54] = _psd_monitor_end_var._position_absolute;
  instrument->_position_relative[54] = _psd_monitor_end_var._position_relative;
    _psd_monitor_end_var._position_relative_is_zero =  coords_test_zero(_psd_monitor_end_var._position_relative);
  instrument->counter_N[54]  = instrument->counter_P[54] = instrument->counter_P2[54] = 0;
  instrument->counter_AbsorbProp[54]= 0;
  return(0);
} /* _psd_monitor_end_setpos */

/* component E_PSD_mon_end=Monitor_nD() SETTING, POSITION/ROTATION */
int _E_PSD_mon_end_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_E_PSD_mon_end_setpos] component E_PSD_mon_end=Monitor_nD() SETTING [C:\\mcstas-3.4\\lib\\monitors\\Monitor_nD.comp:261]");
  stracpy(_E_PSD_mon_end_var._name, "E_PSD_mon_end", 16384);
  stracpy(_E_PSD_mon_end_var._type, "Monitor_nD", 16384);
  _E_PSD_mon_end_var._index=55;
  int current_setpos_index = 55;
  if("" && strlen(""))
    stracpy(_E_PSD_mon_end_var._parameters.user1, "" ? "" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.user1[0]='\0';
  if("" && strlen(""))
    stracpy(_E_PSD_mon_end_var._parameters.user2, "" ? "" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.user2[0]='\0';
  if("" && strlen(""))
    stracpy(_E_PSD_mon_end_var._parameters.user3, "" ? "" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.user3[0]='\0';
  _E_PSD_mon_end_var._parameters.xwidth = _instrument_var._parameters.Ld;
  _E_PSD_mon_end_var._parameters.yheight = 0.5;
  _E_PSD_mon_end_var._parameters.zdepth = 0;
  _E_PSD_mon_end_var._parameters.xmin = 0;
  _E_PSD_mon_end_var._parameters.xmax = 0;
  _E_PSD_mon_end_var._parameters.ymin = 0;
  _E_PSD_mon_end_var._parameters.ymax = 0;
  _E_PSD_mon_end_var._parameters.zmin = 0;
  _E_PSD_mon_end_var._parameters.zmax = 0;
  _E_PSD_mon_end_var._parameters.bins = 0;
  _E_PSD_mon_end_var._parameters.min = -1e40;
  _E_PSD_mon_end_var._parameters.max = 1e40;
  _E_PSD_mon_end_var._parameters.restore_neutron = 1;
  _E_PSD_mon_end_var._parameters.radius = 0;
  if("x bins=1001 energy limits=[1 7] bins=1201" && strlen("x bins=1001 energy limits=[1 7] bins=1201"))
    stracpy(_E_PSD_mon_end_var._parameters.options, "x bins=1001 energy limits=[1 7] bins=1201" ? "x bins=1001 energy limits=[1 7] bins=1201" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.options[0]='\0';
  if("E_end" && strlen("E_end"))
    stracpy(_E_PSD_mon_end_var._parameters.filename, "E_end" ? "E_end" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.filename[0]='\0';
  if("NULL" && strlen("NULL"))
    stracpy(_E_PSD_mon_end_var._parameters.geometry, "NULL" ? "NULL" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.geometry[0]='\0';
  _E_PSD_mon_end_var._parameters.nowritefile = 0;
  if("NULL" && strlen("NULL"))
    stracpy(_E_PSD_mon_end_var._parameters.username1, "NULL" ? "NULL" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.username1[0]='\0';
  if("NULL" && strlen("NULL"))
    stracpy(_E_PSD_mon_end_var._parameters.username2, "NULL" ? "NULL" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.username2[0]='\0';
  if("NULL" && strlen("NULL"))
    stracpy(_E_PSD_mon_end_var._parameters.username3, "NULL" ? "NULL" : "", 16384);
  else 
  _E_PSD_mon_end_var._parameters.username3[0]='\0';


  /* component E_PSD_mon_end=Monitor_nD() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _detector_var._rotation_absolute, _E_PSD_mon_end_var._rotation_absolute);
    rot_transpose(_psd_monitor_end_var._rotation_absolute, tr1);
    rot_mul(_E_PSD_mon_end_var._rotation_absolute, tr1, _E_PSD_mon_end_var._rotation_relative);
    _E_PSD_mon_end_var._rotation_is_identity =  rot_test_identity(_E_PSD_mon_end_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_detector_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _E_PSD_mon_end_var._position_absolute = coords_add(_detector_var._position_absolute, tc2);
    tc1 = coords_sub(_psd_monitor_end_var._position_absolute, _E_PSD_mon_end_var._position_absolute);
    _E_PSD_mon_end_var._position_relative = rot_apply(_E_PSD_mon_end_var._rotation_absolute, tc1);
  } /* E_PSD_mon_end=Monitor_nD() AT ROTATED */
  DEBUG_COMPONENT("E_PSD_mon_end", _E_PSD_mon_end_var._position_absolute, _E_PSD_mon_end_var._rotation_absolute);
  instrument->_position_absolute[55] = _E_PSD_mon_end_var._position_absolute;
  instrument->_position_relative[55] = _E_PSD_mon_end_var._position_relative;
    _E_PSD_mon_end_var._position_relative_is_zero =  coords_test_zero(_E_PSD_mon_end_var._position_relative);
  instrument->counter_N[55]  = instrument->counter_P[55] = instrument->counter_P2[55] = 0;
  instrument->counter_AbsorbProp[55]= 0;
  return(0);
} /* _E_PSD_mon_end_setpos */

_class_Progress_bar *class_Progress_bar_init(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_init] component origin=Progress_bar() INITIALISE [C:\\mcstas-3.4\\lib\\misc\\Progress_bar.comp:57]");

IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", instrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
  #ifdef OPENACC
  time(&StartTime);
  #endif
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_init */

_class_Source_div *class_Source_div_init(_class_Source_div *_comp
) {
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define gauss (_comp->_parameters.gauss)
  #define flux (_comp->_parameters.flux)
  #define sigmah (_comp->_parameters.sigmah)
  #define sigmav (_comp->_parameters.sigmav)
  #define p_init (_comp->_parameters.p_init)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  SIG_MESSAGE("[_source_init] component source=Source_div() INITIALISE [C:\\mcstas-3.4\\lib\\sources\\Source_div.comp:77]");

sigmah = DEG2RAD*focus_aw/(sqrt(8.0*log(2.0)));
  sigmav = DEG2RAD*focus_ah/(sqrt(8.0*log(2.0)));

  if (xwidth < 0 || yheight < 0 || focus_aw < 0 || focus_ah < 0) {
      printf("Source_div: %s: Error in input parameter values!\n"
             "ERROR       Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }
  if ((!lambda0 && !E0 && !dE && !dlambda)) {
    printf("Source_div: %s: You must specify either a wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(0);
  }
  if ((!lambda0 && !dlambda && (E0 <= 0 || dE < 0 || E0-dE <= 0))
    || (!E0 && !dE && (lambda0 <= 0 || dlambda < 0 || lambda0-dlambda <= 0))) {
    printf("Source_div: %s: Unmeaningful definition of wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }
  /* compute distance to next component */
  Coords ToTarget;
  double tx,ty,tz;
  ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+1),POS_A_CURRENT_COMP);
  ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
  coords_get(ToTarget, &tx, &ty, &tz);
  dist=sqrt(tx*tx+ty*ty+tz*tz);
  /* compute target area */
  if (dist) {
    focus_xw=dist*tan(focus_aw*DEG2RAD);
    focus_yh=dist*tan(focus_ah*DEG2RAD);
  }

  p_init  = flux*1e4*xwidth*yheight/mcget_ncount();
  if (!focus_aw || !focus_ah)
    exit(printf("Source_div: %s: Zero divergence defined. \n"
                "ERROR       Use non zero values for focus_aw and focus_ah.\n",
           NAME_CURRENT_COMP));
  p_init *= 2*fabs(DEG2RAD*focus_aw*sin(DEG2RAD*focus_ah/2));  /* solid angle */
  if (dlambda)
    p_init *= 2*dlambda;
  else if (dE)
    p_init *= 2*dE;
  #undef xwidth
  #undef yheight
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef gauss
  #undef flux
  #undef sigmah
  #undef sigmav
  #undef p_init
  #undef dist
  #undef focus_xw
  #undef focus_yh
  return(_comp);
} /* class_Source_div_init */

_class_PSD_monitor *class_PSD_monitor_init(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_psd_test_init] component psd_test=PSD_monitor() INITIALISE [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:62]");

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(0);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  // Use instance name for monitor output if no input was given
  if (!strcmp(filename,"\0")) sprintf(filename,NAME_CURRENT_COMP);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_init */

_class_E_monitor *class_E_monitor_init(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_init] component e_monitor=E_monitor() INITIALISE [C:\\mcstas-3.4\\lib\\monitors\\E_monitor.comp:69]");

  int i;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("E_monitor: %s: Null detection area !\n"
           "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  E_N = create_darr1d(nE);
  E_p = create_darr1d(nE);
  E_p2 = create_darr1d(nE);

  for (i=0; i<nE; i++)
  {
    E_N[i] = 0;
    E_p[i] = 0;
    E_p2[i] = 0;
  }
  S_p = S_pE = S_pE2 = 0;

  // Use instance name for monitor output if no input was given
  if (!strcmp(filename,"\0")) sprintf(filename,NAME_CURRENT_COMP);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_init */

_class_Monochromator_bent *class_Monochromator_bent_init(_class_Monochromator_bent *_comp
) {
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define xthickness (_comp->_parameters.xthickness)
  #define radius_x (_comp->_parameters.radius_x)
  #define plane_of_reflection (_comp->_parameters.plane_of_reflection)
  #define angle_to_cut_horizontal (_comp->_parameters.angle_to_cut_horizontal)
  #define angle_to_cut_vertical (_comp->_parameters.angle_to_cut_vertical)
  #define mosaicity (_comp->_parameters.mosaicity)
  #define domainthickness (_comp->_parameters.domainthickness)
  #define temperature (_comp->_parameters.temperature)
  #define verbose (_comp->_parameters.verbose)
  #define neutron (_comp->_parameters.neutron)
  #define monochromator (_comp->_parameters.monochromator)
  #define non_scattered (_comp->_parameters.non_scattered)
  #define scattered (_comp->_parameters.scattered)
  #define non_hit (_comp->_parameters.non_hit)
  SIG_MESSAGE("[_Analyzer_0_init] component Analyzer_0=Monochromator_bent() INITIALISE [Monochromator_bent.comp:739]");

	non_scattered = 0;
	scattered = 0;
	non_hit = 0;
/////////////////////
// ERROR FUNCTIONS //
/////////////////////
	if (radius_x < 0)
		exit(printf("G_approach_bent_perfect_crystal: %s: incorrect radius_x=%g\n", NAME_CURRENT_COMP, radius_x));
	if (xthickness <= 0)
		exit(printf("G_approach_bent_perfect_crystal: %s: invalid monochromator xthickness=%g\n", NAME_CURRENT_COMP, xthickness));
	if (zwidth <= 0)
		exit(printf("G_approach_bent_perfect_crystal: %s: invalid monochromator zwidth=%g\n", NAME_CURRENT_COMP, zwidth));

//////////////////////////
// Monochromator values //
//////////////////////////

	// Which type of monochromator
		if (!radius_x && !mosaicity) {
			monochromator.type = 0;
		} else if (radius_x && !mosaicity) {
			monochromator.type = 1;
		} else if (!radius_x && mosaicity) {
			monochromator.type = 2;
		} else if (radius_x && mosaicity) {
			monochromator.type = 3;
		}
	
		// Define monochromator angles
	if ((monochromator.type == 0) || (monochromator.type == 2)){
		monochromator.max_angle = PI;
		monochromator.min_angle = PI;
	} else {
		monochromator.max_angle = atan(zwidth/(2*radius_x)) + PI;
		monochromator.min_angle = -atan(zwidth/(2*radius_x)) + PI;
	}
	
	// Read the designated plane of reflection, for use in the monochromator
		enum crystal_plane plane = stringToEnum(&plane_of_reflection);
	
	// Set monochromator values
		monochromator.length = zwidth;
		monochromator.height = yheight;
		monochromator.thickness = xthickness;
		monochromator.radius_horizontal = radius_x;
		monochromator.domain_thickness = domainthickness;
		monochromator.lattice_spacing = crystal_table[plane][0];
		monochromator.Maier_Leibnitz_reflectivity = crystal_table[plane][1]*100; // Convert to SI and Angstrom
		monochromator.bound_atom_scattering_cross_section = crystal_table[plane][2];
		monochromator.absorption_for_1AA_Neutrons = crystal_table[plane][3];
		monochromator.incoherent_scattering_cross_section = crystal_table[plane][4];
		monochromator.volume = crystal_table[plane][5];
		monochromator.atomic_number = crystal_table[plane][6];
		monochromator.debye_temperature = crystal_table[plane][7];
		monochromator.Constant_from_Freund_paper = crystal_table[plane][8];
		monochromator.poisson_ratio = crystal_table[plane][9];
		monochromator.temperature_mono = 300;

	// Calculate Debye Waller Factor
		calculate_B0_and_BT(&monochromator);
		monochromator.Debye_Waller_factor = exp(-(monochromator.B0 + monochromator.BT)/2/square(monochromator.lattice_spacing));

	// Set mosaicity
		if ((monochromator.type == 2) || (monochromator.type == 3)){
			//Input mosaicity is in arc min. Convert to radians, then to sigma from FWHM
			monochromator.mos = mosaicity*MIN2RAD*FWHM2RMS; 
		}
		else {
			monochromator.mos = 0; 
		}   

	// Set scattering vector
		angle_to_cut_horizontal *= DEG2RAD;
		angle_to_cut_vertical *= DEG2RAD;

		monochromator.G_size_zero = 2*PI/monochromator.lattice_spacing;
		
		monochromator.G0[0] = monochromator.G_size_zero*cos(angle_to_cut_horizontal)*cos(angle_to_cut_vertical);
		monochromator.G0[1] = monochromator.G_size_zero*sin(angle_to_cut_vertical);
		monochromator.G0[2] = monochromator.G_size_zero*sin(angle_to_cut_horizontal)*cos(angle_to_cut_vertical);

		monochromator.perp_G[0] = sin(angle_to_cut_horizontal)*cos(angle_to_cut_vertical);
		monochromator.perp_G[1] = sin(angle_to_cut_vertical);
		monochromator.perp_G[2] = -cos(angle_to_cut_horizontal)*cos(angle_to_cut_vertical);

	// Initialize lattice_spacing_gradient_field 
		if ((monochromator.type == 0 )|| (monochromator.type == 2)){
			monochromator.lattice_spacing_gradient_field[0][0] = 0;
			monochromator.lattice_spacing_gradient_field[0][1] = 0;
			monochromator.lattice_spacing_gradient_field[0][2] = 0;
			monochromator.lattice_spacing_gradient_field[1][0] = 0;
			monochromator.lattice_spacing_gradient_field[1][1] = 0;
			monochromator.lattice_spacing_gradient_field[1][2] = 0;
			monochromator.lattice_spacing_gradient_field[2][0] = 0;
			monochromator.lattice_spacing_gradient_field[2][1] = 0;
			monochromator.lattice_spacing_gradient_field[2][2] = 0;
		} else{
			double curvature = 1/radius_x;
			monochromator.lattice_spacing_gradient_field[0][0] = -monochromator.poisson_ratio*cos(angle_to_cut_horizontal)*monochromator.G_size_zero*curvature;
			monochromator.lattice_spacing_gradient_field[0][1] = 0;
			monochromator.lattice_spacing_gradient_field[0][2] = sin(angle_to_cut_horizontal)*monochromator.G_size_zero*curvature;
			monochromator.lattice_spacing_gradient_field[1][0] = 0;
			monochromator.lattice_spacing_gradient_field[1][1] = 0;
			monochromator.lattice_spacing_gradient_field[1][2] = 0;
			monochromator.lattice_spacing_gradient_field[2][0] = sin(angle_to_cut_horizontal)*monochromator.G_size_zero*curvature;
			monochromator.lattice_spacing_gradient_field[2][1] = 0;
			monochromator.lattice_spacing_gradient_field[2][2] = -cos(angle_to_cut_horizontal)*monochromator.G_size_zero*curvature;
		}
  #undef zwidth
  #undef yheight
  #undef xthickness
  #undef radius_x
  #undef plane_of_reflection
  #undef angle_to_cut_horizontal
  #undef angle_to_cut_vertical
  #undef mosaicity
  #undef domainthickness
  #undef temperature
  #undef verbose
  #undef neutron
  #undef monochromator
  #undef non_scattered
  #undef scattered
  #undef non_hit
  return(_comp);
} /* class_Monochromator_bent_init */

_class_Monitor_nD *class_Monitor_nD_init(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  SIG_MESSAGE("[_E_PSD_mon_end_init] component E_PSD_mon_end=Monitor_nD() INITIALISE [C:\\mcstas-3.4\\lib\\monitors\\Monitor_nD.comp:261]");

  char tmp[CHAR_BUF_LENGTH];
  strcpy(Vars.compcurname, NAME_CURRENT_COMP);
  if (options != NULL)
    strncpy(Vars.option, options, CHAR_BUF_LENGTH);
  else {
    strcpy(Vars.option, "x y");
    printf("Monitor_nD: %s has no option specified. Setting to PSD ('x y') monitor.\n", NAME_CURRENT_COMP);
  }
  Vars.compcurpos = POS_A_CURRENT_COMP;

  if (strstr(Vars.option, "source"))
    strcat(Vars.option, " list, x y z vx vy vz t sx sy sz ");

  if (bins) { sprintf(tmp, " all bins=%ld ", (long)bins); strcat(Vars.option, tmp); }
  if (min > -FLT_MAX && max < FLT_MAX) { sprintf(tmp, " all limits=[%g %g]", min, max); strcat(Vars.option, tmp); }
  else if (min > -FLT_MAX) { sprintf(tmp, " all min=%g", min); strcat(Vars.option, tmp); }
  else if (max <  FLT_MAX) { sprintf(tmp, " all max=%g", max); strcat(Vars.option, tmp); }

  /* transfer, "zero", and check username- and user variable strings to Vars struct*/
  strncpy(Vars.UserName1,
    username1 && strlen(username1) && strcmp(username1, "0") && strcmp(username1, "NULL") ?
    username1 : "", 128);
  strncpy(Vars.UserName2,
    username2 && strlen(username2) && strcmp(username2, "0") && strcmp(username2, "NULL") ?
    username2 : "", 128);
  strncpy(Vars.UserName3,
    username3 && strlen(username3) && strcmp(username3, "0") && strcmp(username3, "NULL") ?
    username3 : "", 128);
  if(user1 && strlen(user1) && strcmp(user1, "0") && strcmp(user1, "NULL")){
    strncpy(Vars.UserVariable1,user1,128);
    int fail;_class_particle testparticle;
    particle_getvar(&testparticle,Vars.UserVariable1,&fail);
    if(fail){
      fprintf(stderr,"Warning (%s): user1=%s is unknown. The signal will not be resolved - this is likely not what you intended.\n",NAME_CURRENT_COMP,user1);
    }
  }
  if(user2 && strlen(user2) && strcmp(user2, "0") && strcmp(user2, "NULL")){
    strncpy(Vars.UserVariable2,user2,128);
    int fail;_class_particle testparticle;
    particle_getvar(&testparticle,Vars.UserVariable2,&fail);
    if(fail){
      fprintf(stderr,"Warning (%s): user2=%s is unknown. The signal will not be resolved - this is likely not what you intended.\n",NAME_CURRENT_COMP,user2);
    }
  }
  if(user3 && strlen(user3) && strcmp(user3, "0") && strcmp(user3, "NULL")){
    strncpy(Vars.UserVariable3,user3,128);
    int fail;_class_particle testparticle;
    particle_getvar(&testparticle,Vars.UserVariable3,&fail);
    if(fail){
      fprintf(stderr,"Warning (%s): user3=%s is unknown. The signal will not be resolved - this is likely not what you intended.\n",NAME_CURRENT_COMP,user3);
    }
  }
 
  /*sanitize parameters set for curved shapes*/
  if(strstr(Vars.option,"cylinder") || strstr(Vars.option,"banana") || strstr(Vars.option,"sphere")){
    /*this _is_ an explicit curved shape. Should have a radius. Inherit from xwidth or zdepth (diameters), x has precedence.*/
    if (!radius){
      if(xwidth){
	radius=xwidth/2.0;
      }else{
	radius=zdepth/2.0;
      }
    }else{
      xwidth=2*radius;
    }
    if(!yheight){
      /*if not set - use the diameter as height for the curved object. This will likely only happen for spheres*/
      yheight=2*radius;
    }
  }else if (radius) {
    /*radius is set - this must be a curved shape. Infer shape from yheight, and set remaining values
     (xwidth etc. They are used inside monitor_nd-lib.*/
    xwidth = zdepth = 2*radius;
    if (yheight){
      /*a height is given (and no shape explitly set - assume cylinder*/
      strcat(Vars.option, " banana");
    }else {
      strcat(Vars.option, " sphere");
      yheight=2*radius;
    }
  }

  int offflag=0;
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL")) {
    #ifndef USE_OFF
    fprintf(stderr,"Error: You are attempting to use an OFF geometry without -DUSE_OFF. You will need to recompile with that define set!\n");
    exit(-1);
    #else
    if (!off_init(  geometry, xwidth, yheight, zdepth, 1, &offdata )) {
      printf("Monitor_nD: %s could not initiate the OFF geometry %s. \n"
             "            Defaulting to normal Monitor dimensions.\n",
             NAME_CURRENT_COMP, geometry);
      strcpy(geometry, "");
    } else {
      offflag=1;
    }
    #endif
  }

  if (!radius && !xwidth && !yheight && !zdepth && !xmin && !xmax && !ymin && !ymax &&
    !strstr(Vars.option, "previous") && (!geometry || !strlen(geometry)))
    exit(printf("Monitor_nD: %s has no dimension specified. Aborting (radius, xwidth, yheight, zdepth, previous, geometry).\n", NAME_CURRENT_COMP));

  Monitor_nD_Init(&DEFS, &Vars, xwidth, yheight, zdepth, xmin,xmax,ymin,ymax,zmin,zmax,offflag);

  if (Vars.Flag_OFF) {
    offdata.mantidflag=Vars.Flag_mantid;
    offdata.mantidoffset=Vars.Coord_Min[Vars.Coord_Number-1];
  }


  if (filename && strlen(filename) && strcmp(filename,"NULL") && strcmp(filename,"0"))
    strncpy(Vars.Mon_File, filename, 128);

  /* check if user given filename with ext will be used more than once */
  if ( ((Vars.Flag_Multiple && Vars.Coord_Number > 1) || Vars.Flag_List) && strchr(Vars.Mon_File,'.') )
  { char *XY; XY = strrchr(Vars.Mon_File,'.'); *XY='_'; }

  if (restore_neutron) Vars.Flag_parallel=1;
  detector.m = 0;

#ifdef USE_MPI
MPI_MASTER(
  if (strstr(Vars.option, "auto") && mpi_node_count > 1)
    printf("Monitor_nD: %s is using automatic limits option 'auto' together with MPI.\n"
           "WARNING     this may create incorrect distributions (but integrated flux will be right).\n", NAME_CURRENT_COMP);
);
#else
#ifdef OPENACC
  if (strstr(Vars.option, "auto"))
    printf("Monitor_nD: %s is requesting automatic limits option 'auto' together with OpenACC.\n"
           "WARNING     this feature is NOT supported using OpenACC and has been disabled!\n", NAME_CURRENT_COMP);
#endif
#endif

  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef nowritefile
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_init */



int init(void) { /* called by mccode_main for template_simple:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "template_simple", 256);

  _origin_setpos(); /* type Progress_bar */
  _source_setpos(); /* type Source_div */
  _psd_test_setpos(); /* type PSD_monitor */
  _e_monitor_setpos(); /* type E_monitor */
  _Analyzer_0_setpos(); /* type Monochromator_bent */
  _Analyzer_1_setpos(); /* type Monochromator_bent */
  _Analyzer_2_setpos(); /* type Monochromator_bent */
  _Analyzer_3_setpos(); /* type Monochromator_bent */
  _Analyzer_4_setpos(); /* type Monochromator_bent */
  _Analyzer_5_setpos(); /* type Monochromator_bent */
  _Analyzer_6_setpos(); /* type Monochromator_bent */
  _Analyzer_7_setpos(); /* type Monochromator_bent */
  _Analyzer_8_setpos(); /* type Monochromator_bent */
  _Analyzer_9_setpos(); /* type Monochromator_bent */
  _Analyzer_10_setpos(); /* type Monochromator_bent */
  _Analyzer_11_setpos(); /* type Monochromator_bent */
  _Analyzer_12_setpos(); /* type Monochromator_bent */
  _Analyzer_13_setpos(); /* type Monochromator_bent */
  _Analyzer_14_setpos(); /* type Monochromator_bent */
  _Analyzer_15_setpos(); /* type Monochromator_bent */
  _Analyzer_16_setpos(); /* type Monochromator_bent */
  _Analyzer_17_setpos(); /* type Monochromator_bent */
  _Analyzer_18_setpos(); /* type Monochromator_bent */
  _Analyzer_19_setpos(); /* type Monochromator_bent */
  _Analyzer_20_setpos(); /* type Monochromator_bent */
  _Analyzer_21_setpos(); /* type Monochromator_bent */
  _Analyzer_22_setpos(); /* type Monochromator_bent */
  _Analyzer_23_setpos(); /* type Monochromator_bent */
  _Analyzer_24_setpos(); /* type Monochromator_bent */
  _Analyzer_25_setpos(); /* type Monochromator_bent */
  _Analyzer_26_setpos(); /* type Monochromator_bent */
  _Analyzer_27_setpos(); /* type Monochromator_bent */
  _Analyzer_28_setpos(); /* type Monochromator_bent */
  _Analyzer_29_setpos(); /* type Monochromator_bent */
  _Analyzer_30_setpos(); /* type Monochromator_bent */
  _Analyzer_31_setpos(); /* type Monochromator_bent */
  _Analyzer_32_setpos(); /* type Monochromator_bent */
  _Analyzer_33_setpos(); /* type Monochromator_bent */
  _Analyzer_34_setpos(); /* type Monochromator_bent */
  _Analyzer_35_setpos(); /* type Monochromator_bent */
  _Analyzer_36_setpos(); /* type Monochromator_bent */
  _Analyzer_37_setpos(); /* type Monochromator_bent */
  _Analyzer_38_setpos(); /* type Monochromator_bent */
  _Analyzer_39_setpos(); /* type Monochromator_bent */
  _Analyzer_40_setpos(); /* type Monochromator_bent */
  _Analyzer_41_setpos(); /* type Monochromator_bent */
  _Analyzer_42_setpos(); /* type Monochromator_bent */
  _Analyzer_43_setpos(); /* type Monochromator_bent */
  _Analyzer_44_setpos(); /* type Monochromator_bent */
  _Analyzer_45_setpos(); /* type Monochromator_bent */
  _rotator_setpos(); /* type Arm */
  _detector_pos_setpos(); /* type Arm */
  _detector_setpos(); /* type Arm */
  _psd_monitor_end_setpos(); /* type PSD_monitor */
  _E_PSD_mon_end_setpos(); /* type Monitor_nD */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(&_origin_var);

  class_Source_div_init(&_source_var);

  class_PSD_monitor_init(&_psd_test_var);

  class_E_monitor_init(&_e_monitor_var);

  class_Monochromator_bent_init(&_Analyzer_0_var);

  class_Monochromator_bent_init(&_Analyzer_1_var);

  class_Monochromator_bent_init(&_Analyzer_2_var);

  class_Monochromator_bent_init(&_Analyzer_3_var);

  class_Monochromator_bent_init(&_Analyzer_4_var);

  class_Monochromator_bent_init(&_Analyzer_5_var);

  class_Monochromator_bent_init(&_Analyzer_6_var);

  class_Monochromator_bent_init(&_Analyzer_7_var);

  class_Monochromator_bent_init(&_Analyzer_8_var);

  class_Monochromator_bent_init(&_Analyzer_9_var);

  class_Monochromator_bent_init(&_Analyzer_10_var);

  class_Monochromator_bent_init(&_Analyzer_11_var);

  class_Monochromator_bent_init(&_Analyzer_12_var);

  class_Monochromator_bent_init(&_Analyzer_13_var);

  class_Monochromator_bent_init(&_Analyzer_14_var);

  class_Monochromator_bent_init(&_Analyzer_15_var);

  class_Monochromator_bent_init(&_Analyzer_16_var);

  class_Monochromator_bent_init(&_Analyzer_17_var);

  class_Monochromator_bent_init(&_Analyzer_18_var);

  class_Monochromator_bent_init(&_Analyzer_19_var);

  class_Monochromator_bent_init(&_Analyzer_20_var);

  class_Monochromator_bent_init(&_Analyzer_21_var);

  class_Monochromator_bent_init(&_Analyzer_22_var);

  class_Monochromator_bent_init(&_Analyzer_23_var);

  class_Monochromator_bent_init(&_Analyzer_24_var);

  class_Monochromator_bent_init(&_Analyzer_25_var);

  class_Monochromator_bent_init(&_Analyzer_26_var);

  class_Monochromator_bent_init(&_Analyzer_27_var);

  class_Monochromator_bent_init(&_Analyzer_28_var);

  class_Monochromator_bent_init(&_Analyzer_29_var);

  class_Monochromator_bent_init(&_Analyzer_30_var);

  class_Monochromator_bent_init(&_Analyzer_31_var);

  class_Monochromator_bent_init(&_Analyzer_32_var);

  class_Monochromator_bent_init(&_Analyzer_33_var);

  class_Monochromator_bent_init(&_Analyzer_34_var);

  class_Monochromator_bent_init(&_Analyzer_35_var);

  class_Monochromator_bent_init(&_Analyzer_36_var);

  class_Monochromator_bent_init(&_Analyzer_37_var);

  class_Monochromator_bent_init(&_Analyzer_38_var);

  class_Monochromator_bent_init(&_Analyzer_39_var);

  class_Monochromator_bent_init(&_Analyzer_40_var);

  class_Monochromator_bent_init(&_Analyzer_41_var);

  class_Monochromator_bent_init(&_Analyzer_42_var);

  class_Monochromator_bent_init(&_Analyzer_43_var);

  class_Monochromator_bent_init(&_Analyzer_44_var);

  class_Monochromator_bent_init(&_Analyzer_45_var);




  class_PSD_monitor_init(&_psd_monitor_end_var);

  class_Monitor_nD_init(&_E_PSD_mon_end_var);

  if (mcdotrace) display();
  DEBUG_INSTR_END();

#ifdef OPENACC
#include <openacc.h>
#pragma acc update device(_origin_var)
#pragma acc update device(_source_var)
#pragma acc update device(_psd_test_var)
#pragma acc update device(_e_monitor_var)
#pragma acc update device(_Analyzer_0_var)
#pragma acc update device(_Analyzer_1_var)
#pragma acc update device(_Analyzer_2_var)
#pragma acc update device(_Analyzer_3_var)
#pragma acc update device(_Analyzer_4_var)
#pragma acc update device(_Analyzer_5_var)
#pragma acc update device(_Analyzer_6_var)
#pragma acc update device(_Analyzer_7_var)
#pragma acc update device(_Analyzer_8_var)
#pragma acc update device(_Analyzer_9_var)
#pragma acc update device(_Analyzer_10_var)
#pragma acc update device(_Analyzer_11_var)
#pragma acc update device(_Analyzer_12_var)
#pragma acc update device(_Analyzer_13_var)
#pragma acc update device(_Analyzer_14_var)
#pragma acc update device(_Analyzer_15_var)
#pragma acc update device(_Analyzer_16_var)
#pragma acc update device(_Analyzer_17_var)
#pragma acc update device(_Analyzer_18_var)
#pragma acc update device(_Analyzer_19_var)
#pragma acc update device(_Analyzer_20_var)
#pragma acc update device(_Analyzer_21_var)
#pragma acc update device(_Analyzer_22_var)
#pragma acc update device(_Analyzer_23_var)
#pragma acc update device(_Analyzer_24_var)
#pragma acc update device(_Analyzer_25_var)
#pragma acc update device(_Analyzer_26_var)
#pragma acc update device(_Analyzer_27_var)
#pragma acc update device(_Analyzer_28_var)
#pragma acc update device(_Analyzer_29_var)
#pragma acc update device(_Analyzer_30_var)
#pragma acc update device(_Analyzer_31_var)
#pragma acc update device(_Analyzer_32_var)
#pragma acc update device(_Analyzer_33_var)
#pragma acc update device(_Analyzer_34_var)
#pragma acc update device(_Analyzer_35_var)
#pragma acc update device(_Analyzer_36_var)
#pragma acc update device(_Analyzer_37_var)
#pragma acc update device(_Analyzer_38_var)
#pragma acc update device(_Analyzer_39_var)
#pragma acc update device(_Analyzer_40_var)
#pragma acc update device(_Analyzer_41_var)
#pragma acc update device(_Analyzer_42_var)
#pragma acc update device(_Analyzer_43_var)
#pragma acc update device(_Analyzer_44_var)
#pragma acc update device(_Analyzer_45_var)
#pragma acc update device(_rotator_var)
#pragma acc update device(_detector_pos_var)
#pragma acc update device(_detector_var)
#pragma acc update device(_psd_monitor_end_var)
#pragma acc update device(_E_PSD_mon_end_var)
#pragma acc update device(_instrument_var)
#endif

  return(0);
} /* init */

/*******************************************************************************
* components TRACE
*******************************************************************************/

#define x (_particle->x)
#define y (_particle->y)
#define z (_particle->z)
#define vx (_particle->vx)
#define vy (_particle->vy)
#define vz (_particle->vz)
#define t (_particle->t)
#define sx (_particle->sx)
#define sy (_particle->sy)
#define sz (_particle->sz)
#define p (_particle->p)
#define mcgravitation (_particle->mcgravitation)
#define mcMagnet (_particle->mcMagnet)
#define allow_backprop (_particle->allow_backprop)
#define _mctmp_a (_particle->_mctmp_a)
#define _mctmp_b (_particle->_mctmp_b)
#define _mctmp_c (_particle->_mctmp_c)
/* if on GPU, globally nullify sprintf,fprintf,printfs   */
/* (Similar defines are available in each comp trace but */
/*  those are not enough to handle external libs etc. )  */
#ifdef OPENACC
#ifndef MULTICORE
#define fprintf(stderr,...) printf(__VA_ARGS__)
#define sprintf(string,...) printf(__VA_ARGS__)
#define exit(...) noprintf()
#define strcmp(a,b) str_comp(a,b)
#define strlen(a) str_len(a)
#endif
#endif
#define SCATTERED (_particle->_scattered)
#define RESTORE (_particle->_restore)
#define RESTORE_NEUTRON(_index, ...) _particle->_restore = _index;
#define ABSORBED (_particle->_absorbed)
#define mcget_run_num() _particle->_uid
#define ABSORB0 do { DEBUG_STATE(); DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; return(_comp); } while(0)
#define ABSORB ABSORB0
#pragma acc routine
_class_Progress_bar *class_Progress_bar_trace(_class_Progress_bar *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_trace] component origin=Progress_bar() TRACE [C:\\mcstas-3.4\\lib\\misc\\Progress_bar.comp:73]");

#ifndef OPENACC
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10 && ncount) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  /* display percentage when percent or minutes have reached step */
  if (EndTime && mcget_ncount() &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;

    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) save(NULL);
  }
#endif
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_trace */

#pragma acc routine
_class_Source_div *class_Source_div_trace(_class_Source_div *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define gauss (_comp->_parameters.gauss)
  #define flux (_comp->_parameters.flux)
  #define sigmah (_comp->_parameters.sigmah)
  #define sigmav (_comp->_parameters.sigmav)
  #define p_init (_comp->_parameters.p_init)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  SIG_MESSAGE("[_source_trace] component source=Source_div() TRACE [C:\\mcstas-3.4\\lib\\sources\\Source_div.comp:123]");

  double E,lambda,v;
  double tan_h;
  double tan_v;
  double thetah;
  double thetav;

  p=p_init;
  z=0;
  t=0;

  x=randpm1()*xwidth/2.0;
  y=randpm1()*yheight/2.0;
  if(lambda0==0) {
    if (!gauss) {
      E=E0+dE*randpm1();              /*  Choose from uniform distribution */
    } else {
      E=E0+randnorm()*dE;
    }
    v=sqrt(E)*SE2V;
  } else {
    if (!gauss) {
      lambda=lambda0+dlambda*randpm1();
    } else {
      lambda=lambda0+randnorm()*dlambda;
    }
    v = K2V*(2*PI/lambda);
  }

  if (gauss==1) {
    thetah = randnorm()*sigmah;
    thetav = randnorm()*sigmav;
  } else {
    thetah = randpm1()*focus_aw*DEG2RAD/2;
    thetav = randpm1()*focus_ah*DEG2RAD/2;
  }

  tan_h = tan(thetah);
  tan_v = tan(thetav);

  /* Perform the correct treatment - no small angle approx. here! */
  vz = v / sqrt(1 + tan_v*tan_v + tan_h*tan_h);
  vy = tan_v * vz;
  vx = tan_h * vz;
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef xwidth
  #undef yheight
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef gauss
  #undef flux
  #undef sigmah
  #undef sigmav
  #undef p_init
  #undef dist
  #undef focus_xw
  #undef focus_yh
  return(_comp);
} /* class_Source_div_trace */

#pragma acc routine
_class_PSD_monitor *class_PSD_monitor_trace(_class_PSD_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_psd_test_trace] component psd_test=PSD_monitor() TRACE [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:82]");

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));

    double p2 = p*p;
    #pragma acc atomic
    PSD_N[i][j] = PSD_N[i][j]+1;

    #pragma acc atomic
    PSD_p[i][j] = PSD_p[i][j]+p;
    
    #pragma acc atomic
    PSD_p2[i][j] = PSD_p2[i][j] + p2;
    
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_trace */

#pragma acc routine
_class_E_monitor *class_E_monitor_trace(_class_E_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_trace] component e_monitor=E_monitor() TRACE [C:\\mcstas-3.4\\lib\\monitors\\E_monitor.comp:99]");

  int i;
  double E;

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    E = VS2E*(vx*vx + vy*vy + vz*vz);

    S_p += p;
    S_pE += p*E;
    S_pE2 += p*E*E;

    i = floor((E-Emin)*nE/(Emax-Emin));
    if(i >= 0 && i < nE)
    {
      double p2 = p*p;
      #pragma acc atomic
      E_N[i] = E_N[i] +1;
      #pragma acc atomic
      E_p[i] = E_p[i] + p;
      #pragma acc atomic
      E_p2[i] = E_p2[i] + p2;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_trace */

#pragma acc routine
_class_Monochromator_bent *class_Monochromator_bent_trace(_class_Monochromator_bent *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define xthickness (_comp->_parameters.xthickness)
  #define radius_x (_comp->_parameters.radius_x)
  #define plane_of_reflection (_comp->_parameters.plane_of_reflection)
  #define angle_to_cut_horizontal (_comp->_parameters.angle_to_cut_horizontal)
  #define angle_to_cut_vertical (_comp->_parameters.angle_to_cut_vertical)
  #define mosaicity (_comp->_parameters.mosaicity)
  #define domainthickness (_comp->_parameters.domainthickness)
  #define temperature (_comp->_parameters.temperature)
  #define verbose (_comp->_parameters.verbose)
  #define neutron (_comp->_parameters.neutron)
  #define monochromator (_comp->_parameters.monochromator)
  #define non_scattered (_comp->_parameters.non_scattered)
  #define scattered (_comp->_parameters.scattered)
  #define non_hit (_comp->_parameters.non_hit)
  SIG_MESSAGE("[_Analyzer_0_trace] component Analyzer_0=Monochromator_bent() TRACE [Monochromator_bent.comp:851]");

	// Set the neutron values for the new neutron and direction
	set_neutron_values(&neutron,x,y,z,vx,vy,vz);
	neutron.direction = 1;
	neutron.scatter_count = 0;
	// Check weight is not 0
	double weight_init = p;
	double k_initial = neutron.ki_size;
									 
	if (weight_init == 0.0){
		ABSORB;
	}
	
	// Find the intersection 
	calc_intersection(&monochromator,&neutron);
																					 
																						

	// Propogate neutron and reset the neutron values
	PROP_DT(neutron.entry_time);
	set_neutron_values(&neutron,x,y,z,vx,vy,vz);
	double angle_on_inner_cylinder;
	// Is the neutron in the area of the monochromator?
	if ((monochromator.type== 1) || (monochromator.type == 3)){
		angle_on_inner_cylinder = PI - atan(neutron.r[2]/monochromator.radius_horizontal);
	} else {
		angle_on_inner_cylinder = PI;
	}
	if ((monochromator.max_angle >= angle_on_inner_cylinder) && (angle_on_inner_cylinder >= monochromator.min_angle)){
		// If yes, then it can possibly be scattered

		// Set the maximal tau for first scattering
		neutron.max_tau = neutron.exit_time-neutron.entry_time;
		neutron.path = 0;
		// Allow for up to 5 scattering indicents to limit computation time
		for (int i = 0; i < 1; i++){	
			// Walk step is repeated for each step here.
			// First calculate scattering vector, eps_zero and beta with random vertical mosaicity
			//prop_half(&monochromator,&neutron);
			calculate_G(&monochromator,&neutron);
			calculate_epszero_and_beta(&monochromator,&neutron);
			// Calculate probability of scattering and time of scattering
			calculate_tau_and_prob(&monochromator,&neutron,i);
			// Will the neutron scatter inside the crystal (between entry and exit times) 
			// and is the probability not INCREDIBLY low (limit calculations)
			if ((neutron.tau[i] < neutron.max_tau) && (neutron.tau[i] > 0)){ 
				
				// Reflect neutron
				// Propogate to scattering point
				scattered++;
				PROP_DT(neutron.tau[i]);
				set_neutron_values(&neutron,x,y,z,vx,vy,vz);

				// Update probability
				p *= neutron.P[i];
				if (p == 0.0) ABSORB;
				neutron.scatter_count += 1;
				
				

				// Update kf and scatter
				calculate_kf(&monochromator,&neutron,i);
				if (i==0){SCATTER;}

				// Reflect neutron and set new variables for next scattering
				reflect_neutron(&neutron, &vx, &vy, &vz);
				set_neutron_values(&neutron,x,y,z,vx,vy,vz);

				neutron.path += neutron.tau[i]*neutron.v_size;
			}else {
				non_scattered++;
				// Travel through the entire path and do not update any speeds
				p *= 1-neutron.P[i];
				neutron.path += sign(neutron.max_tau)*neutron.max_tau*neutron.v_size;
				// Break out of loop, such that no more scattering can occur
				break;
			}
			// Update direction 
			neutron.direction *= -1;

			// Calculate new intersection points
			calc_intersection(&monochromator,&neutron);
			
			if (neutron.direction == -1){
				neutron.max_tau = neutron.entry_time;
				
			} else {
				neutron.max_tau = neutron.exit_time;
			}
			
		
		}
		// Attenuate probability for path
		double attenuation_coefficient = calculate_attenuation_coefficient(&monochromator, &neutron);
		p *= exp(-attenuation_coefficient*neutron.path);
	}else {
		// If not, let it pass through the component
		/* restore neutron state when no interaction */
		non_hit++;
		RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
	}
	
	// # Done
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif

if (_comp->_index == 5) { // EXTEND 'Analyzer_0'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 6) { // EXTEND 'Analyzer_1'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 7) { // EXTEND 'Analyzer_2'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 8) { // EXTEND 'Analyzer_3'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 9) { // EXTEND 'Analyzer_4'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 10) { // EXTEND 'Analyzer_5'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 11) { // EXTEND 'Analyzer_6'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 12) { // EXTEND 'Analyzer_7'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 13) { // EXTEND 'Analyzer_8'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 14) { // EXTEND 'Analyzer_9'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 15) { // EXTEND 'Analyzer_10'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 16) { // EXTEND 'Analyzer_11'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 17) { // EXTEND 'Analyzer_12'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 18) { // EXTEND 'Analyzer_13'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 19) { // EXTEND 'Analyzer_14'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 20) { // EXTEND 'Analyzer_15'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 21) { // EXTEND 'Analyzer_16'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 22) { // EXTEND 'Analyzer_17'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 23) { // EXTEND 'Analyzer_18'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 24) { // EXTEND 'Analyzer_19'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 25) { // EXTEND 'Analyzer_20'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 26) { // EXTEND 'Analyzer_21'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 27) { // EXTEND 'Analyzer_22'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 28) { // EXTEND 'Analyzer_23'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 29) { // EXTEND 'Analyzer_24'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 30) { // EXTEND 'Analyzer_25'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 31) { // EXTEND 'Analyzer_26'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 32) { // EXTEND 'Analyzer_27'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 33) { // EXTEND 'Analyzer_28'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 34) { // EXTEND 'Analyzer_29'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 35) { // EXTEND 'Analyzer_30'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 36) { // EXTEND 'Analyzer_31'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 37) { // EXTEND 'Analyzer_32'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 38) { // EXTEND 'Analyzer_33'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 39) { // EXTEND 'Analyzer_34'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 40) { // EXTEND 'Analyzer_35'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 41) { // EXTEND 'Analyzer_36'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 42) { // EXTEND 'Analyzer_37'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 43) { // EXTEND 'Analyzer_38'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 44) { // EXTEND 'Analyzer_39'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 45) { // EXTEND 'Analyzer_40'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 46) { // EXTEND 'Analyzer_41'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 47) { // EXTEND 'Analyzer_42'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 48) { // EXTEND 'Analyzer_43'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 49) { // EXTEND 'Analyzer_44'
    if(SCATTERED){scat=1;}
}
if (_comp->_index == 50) { // EXTEND 'Analyzer_45'
    if(SCATTERED){scat=1;}
}

  #undef zwidth
  #undef yheight
  #undef xthickness
  #undef radius_x
  #undef plane_of_reflection
  #undef angle_to_cut_horizontal
  #undef angle_to_cut_vertical
  #undef mosaicity
  #undef domainthickness
  #undef temperature
  #undef verbose
  #undef neutron
  #undef monochromator
  #undef non_scattered
  #undef scattered
  #undef non_hit
  return(_comp);
} /* class_Monochromator_bent_trace */

#pragma acc routine
_class_Monitor_nD *class_Monitor_nD_trace(_class_Monitor_nD *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  SIG_MESSAGE("[_E_PSD_mon_end_trace] component E_PSD_mon_end=Monitor_nD() TRACE [C:\\mcstas-3.4\\lib\\monitors\\Monitor_nD.comp:400]");

  double  transmit_he3=1.0;
  double  multiplier_capture=1.0;
  double  t0 = 0;
  double  t1 = 0;
  int     pp;
  int     intersect   = 0;
  char    Flag_Restore = 0;

  #ifdef OPENACC
  #ifdef USE_OFF
  off_struct thread_offdata = offdata;
  #endif
  #else
  #define thread_offdata offdata
  #endif
  
  /* this is done automatically
    STORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  */
  #ifdef USE_OFF
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    /* determine intersections with object */
    intersect = off_intersect_all(&t0, &t1, NULL, NULL,
				  x,y,z, vx, vy, vz, 0, 0, 0, &thread_offdata );
    if (Vars.Flag_mantid) {
      if(intersect) {
        Vars.OFF_polyidx=thread_offdata.nextintersect;
      } else {
        Vars.OFF_polyidx=-1;
      }
    }
  }
  else
  #endif
    if ( (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE)
            || (abs(Vars.Flag_Shape) == DEFS.SHAPE_DISK) ) /* square xy or disk xy */
  {
    // propagate to xy plane and find intersection
    // make sure the event is recoverable afterwards
    t0 = t;
    ALLOW_BACKPROP;
    PROP_Z0;
    if ( (t>=t0) && (z==0.0) ) // forward propagation to xy plane was successful
    {
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SQUARE)
      {
        // square xy
        intersect = (x>=Vars.mxmin && x<=Vars.mxmax && y>=Vars.mymin && y<=Vars.mymax);
      }
      else
      {
        // disk xy
        intersect = (SQR(x) + SQR(y)) <= SQR(Vars.Sphere_Radius);
      }
    }
    else
    {
      intersect=0;
    }
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) /* sphere */
  {
    intersect = sphere_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius);
  /*      intersect = (intersect && t0 > 0); */
  }
  else if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)) /* cylinder */
  {
    intersect = cylinder_intersect(&t0, &t1, x, y, z, vx, vy, vz, Vars.Sphere_Radius, Vars.Cylinder_Height);
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) /* box */
  {
    intersect = box_intersect(&t0, &t1, x, y, z, vx, vy, vz,
                              fabs(Vars.mxmax-Vars.mxmin), fabs(Vars.mymax-Vars.mymin), fabs(Vars.mzmax-Vars.mzmin));
  }
  else if (abs(Vars.Flag_Shape) == DEFS.SHAPE_PREVIOUS) /* previous comp */
  { intersect = 1; }

  if (intersect)
  {
    if ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND)
     || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA)
     || (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL")) )
    {
      /* check if we have to remove the top/bottom with BANANA shape */
        if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA) {
            if (intersect == 1) { // Entered and left through sides
                if (t0 < 0 && t1 > 0) {
                    t0 = t;  /* neutron was already inside ! */
                }
                if (t1 < 0 && t0 > 0) { /* neutron exit before entering !! */
                    t1 = t;
                }
                /* t0 is now time of incoming intersection with the detection area */
                if ((Vars.Flag_Shape < 0) && (t1 > 0)) {
                    PROP_DT(t1); /* t1 outgoing beam */
                } else {
                    PROP_DT(t0); /* t0 incoming beam */
                }
            } else if (intersect == 3 || intersect == 5) { // Entered from top or bottom, left through side
                if ((Vars.Flag_Shape < 0) && (t1 > 0)) {
                    PROP_DT(t1); /* t1 outgoing beam */
                } else {
                    intersect=0;
                    Flag_Restore=1;
                }
            } else if (intersect == 9 || intersect == 17) { // Entered through side, left from top or bottom
                if ((Vars.Flag_Shape < 0) && (t1 > 0)) {
                    intersect=0;
                    Flag_Restore=1;
                } else {
                    PROP_DT(t0); /* t0 incoming beam */
                }
            } else if (intersect == 13 || intersect == 19) { // Went through top/bottom on entry and exit
                intersect=0;
                Flag_Restore=1;
            } else {
                printf("Cylinder_intersect returned unexpected value %l\n", intersect);
            }
        } else {
            // All other shapes than the BANANA 
            if (t0 < 0 && t1 > 0)
              t0 = t;  /* neutron was already inside ! */
            if (t1 < 0 && t0 > 0) /* neutron exit before entering !! */
              t1 = t;
            /* t0 is now time of incoming intersection with the detection area */
            if ((Vars.Flag_Shape < 0) && (t1 > 0))
              PROP_DT(t1); /* t1 outgoing beam */
            else
              PROP_DT(t0); /* t0 incoming beam */
        }
      
      /* Final test if we are on lid / bottom of banana/sphere */
      if (abs(Vars.Flag_Shape) == DEFS.SHAPE_BANANA || abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) {
        if (Vars.Cylinder_Height && fabs(y) >= Vars.Cylinder_Height/2 - FLT_EPSILON) {
          intersect=0;
          Flag_Restore=1;
        }
      }
    }
  }

  if (intersect)
  {
    /* Now get the data to monitor: current or keep from PreMonitor */
/*    if (Vars.Flag_UsePreMonitor != 1)*/
/*    {*/
/*      Vars.cp  = p;*/
/*      Vars.cx  = x;*/
/*      Vars.cvx = vx;*/
/*      Vars.csx = sx;*/
/*      Vars.cy  = y;*/
/*      Vars.cvy = vy;*/
/*      Vars.csy = sy;*/
/*      Vars.cz  = z;*/
/*      Vars.cvz = vz;*/
/*      Vars.csz = sz;*/
/*      Vars.ct  = t;*/
/*    }*/

    if ((Vars.He3_pressure > 0) && (t1 != t0) && ((abs(Vars.Flag_Shape) == DEFS.SHAPE_SPHERE) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_CYLIND) || (abs(Vars.Flag_Shape) == DEFS.SHAPE_BOX)))
    {
      transmit_he3 = exp(-7.417*Vars.He3_pressure*fabs(t1-t0)*2*PI*K2V);
      /* will monitor the absorbed part */
      p = p * (1-transmit_he3);
    }

    if (Vars.Flag_capture)
    {
      multiplier_capture = V2K*sqrt(vx*vx+vy*vy+vz*vz);
      if (multiplier_capture != 0) multiplier_capture = 2*PI/multiplier_capture; /* lambda. lambda(2200 m/2) = 1.7985 Angs  */
      p = p * multiplier_capture/1.7985;
    }

    pp = Monitor_nD_Trace(&DEFS, &Vars, _particle);
    if (pp==0.0)
    {
      ABSORB;
    }
    else if(pp==1)
    {
      SCATTER;
    }

    /*set weight to undetected part if capture and/or he3_pressure*/
    if (Vars.He3_pressure > 0){
      /* after monitor, only remains 1-p_detect */
      p = p * transmit_he3/(1.0-transmit_he3);
    }

    if (Vars.Flag_capture){
      p = p / multiplier_capture*1.7985;
    }

    if (Vars.Flag_parallel) /* back to neutron state before detection */
      Flag_Restore = 1;
  } /* end if intersection */
  else {
    if (Vars.Flag_Absorb && !Vars.Flag_parallel)
    {
      // restore neutron ray before absorbing for correct mcdisplay
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else Flag_Restore = 1;  /* no intersection, back to previous state */
  }

  if (Flag_Restore)
  {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef nowritefile
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_trace */

/* *****************************************************************************
* instrument 'template_simple' TRACE
***************************************************************************** */

#ifndef FUNNEL
#pragma acc routine
int raytrace(_class_particle* _particle) { /* single event propagation, called by mccode_main for template_simple:TRACE */

  /* init variables and counters for TRACE */
  #undef ABSORB0
  #undef ABSORB
  #define ABSORB0 do { DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; return(ABSORBED);} while(0)
  #define ABSORB ABSORB0
  DEBUG_ENTER();
  DEBUG_STATE();
  _particle->flag_nocoordschange=0; /* Init */
  _class_particle _particle_save;
  /* the main iteration loop for one incoming event */
  while (!ABSORBED) { /* iterate event until absorbed */
    /* send particle event to component instance, one after the other */
    /* begin component origin=Progress_bar() [1] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_origin_var._rotation_is_identity) {
        if(!_origin_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _origin_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_origin_var._position_relative, _origin_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 1) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_origin_var._name);
      DEBUG_STATE();
      class_Progress_bar_trace(&_origin_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component origin [1] */
    /* begin component source=Source_div() [2] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_source_var._rotation_is_identity) {
        if(!_source_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _source_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_source_var._position_relative, _source_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 2) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_source_var._name);
      DEBUG_STATE();
      class_Source_div_trace(&_source_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component source [2] */
    /* begin component psd_test=PSD_monitor() [3] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_psd_test_var._rotation_is_identity) {
        if(!_psd_test_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _psd_test_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_psd_test_var._position_relative, _psd_test_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 3) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_psd_test_var._name);
      DEBUG_STATE();
      class_PSD_monitor_trace(&_psd_test_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component psd_test [3] */
    /* begin component e_monitor=E_monitor() [4] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_e_monitor_var._rotation_is_identity) {
        if(!_e_monitor_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _e_monitor_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_e_monitor_var._position_relative, _e_monitor_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 4) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_e_monitor_var._name);
      DEBUG_STATE();
      class_E_monitor_trace(&_e_monitor_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component e_monitor [4] */
    /* begin component Analyzer_0=Monochromator_bent() [5] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_0_var._rotation_is_identity) {
        if(!_Analyzer_0_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_0_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_0_var._position_relative, _Analyzer_0_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 5) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_Analyzer_0_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_0_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_0 [5] */
    /* begin component Analyzer_1=Monochromator_bent() [6] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_1_var._rotation_is_identity) {
        if(!_Analyzer_1_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_1_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_1_var._position_relative, _Analyzer_1_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 6) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_1_var._position_relative, _Analyzer_1_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_1_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_1_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_1 [6] */
    /* begin component Analyzer_2=Monochromator_bent() [7] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_2_var._rotation_is_identity) {
        if(!_Analyzer_2_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_2_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_2_var._position_relative, _Analyzer_2_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 7) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_2_var._position_relative, _Analyzer_2_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_2_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_2_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_2 [7] */
    /* begin component Analyzer_3=Monochromator_bent() [8] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_3_var._rotation_is_identity) {
        if(!_Analyzer_3_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_3_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_3_var._position_relative, _Analyzer_3_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 8) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_3_var._position_relative, _Analyzer_3_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_3_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_3_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_3 [8] */
    /* begin component Analyzer_4=Monochromator_bent() [9] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_4_var._rotation_is_identity) {
        if(!_Analyzer_4_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_4_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_4_var._position_relative, _Analyzer_4_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 9) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_4_var._position_relative, _Analyzer_4_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_4_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_4_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_4 [9] */
    /* begin component Analyzer_5=Monochromator_bent() [10] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_5_var._rotation_is_identity) {
        if(!_Analyzer_5_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_5_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_5_var._position_relative, _Analyzer_5_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 10) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_5_var._position_relative, _Analyzer_5_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_5_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_5_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_5 [10] */
    /* begin component Analyzer_6=Monochromator_bent() [11] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_6_var._rotation_is_identity) {
        if(!_Analyzer_6_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_6_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_6_var._position_relative, _Analyzer_6_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 11) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_6_var._position_relative, _Analyzer_6_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_6_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_6_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_6 [11] */
    /* begin component Analyzer_7=Monochromator_bent() [12] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_7_var._rotation_is_identity) {
        if(!_Analyzer_7_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_7_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_7_var._position_relative, _Analyzer_7_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 12) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_7_var._position_relative, _Analyzer_7_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_7_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_7_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_7 [12] */
    /* begin component Analyzer_8=Monochromator_bent() [13] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_8_var._rotation_is_identity) {
        if(!_Analyzer_8_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_8_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_8_var._position_relative, _Analyzer_8_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 13) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_8_var._position_relative, _Analyzer_8_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_8_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_8_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_8 [13] */
    /* begin component Analyzer_9=Monochromator_bent() [14] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_9_var._rotation_is_identity) {
        if(!_Analyzer_9_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_9_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_9_var._position_relative, _Analyzer_9_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 14) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_9_var._position_relative, _Analyzer_9_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_9_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_9_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_9 [14] */
    /* begin component Analyzer_10=Monochromator_bent() [15] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_10_var._rotation_is_identity) {
        if(!_Analyzer_10_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_10_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_10_var._position_relative, _Analyzer_10_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 15) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_10_var._position_relative, _Analyzer_10_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_10_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_10_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_10 [15] */
    /* begin component Analyzer_11=Monochromator_bent() [16] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_11_var._rotation_is_identity) {
        if(!_Analyzer_11_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_11_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_11_var._position_relative, _Analyzer_11_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 16) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_11_var._position_relative, _Analyzer_11_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_11_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_11_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_11 [16] */
    /* begin component Analyzer_12=Monochromator_bent() [17] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_12_var._rotation_is_identity) {
        if(!_Analyzer_12_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_12_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_12_var._position_relative, _Analyzer_12_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 17) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_12_var._position_relative, _Analyzer_12_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_12_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_12_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_12 [17] */
    /* begin component Analyzer_13=Monochromator_bent() [18] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_13_var._rotation_is_identity) {
        if(!_Analyzer_13_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_13_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_13_var._position_relative, _Analyzer_13_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 18) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_13_var._position_relative, _Analyzer_13_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_13_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_13_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_13 [18] */
    /* begin component Analyzer_14=Monochromator_bent() [19] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_14_var._rotation_is_identity) {
        if(!_Analyzer_14_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_14_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_14_var._position_relative, _Analyzer_14_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 19) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_14_var._position_relative, _Analyzer_14_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_14_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_14_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_14 [19] */
    /* begin component Analyzer_15=Monochromator_bent() [20] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_15_var._rotation_is_identity) {
        if(!_Analyzer_15_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_15_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_15_var._position_relative, _Analyzer_15_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 20) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_15_var._position_relative, _Analyzer_15_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_15_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_15_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_15 [20] */
    /* begin component Analyzer_16=Monochromator_bent() [21] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_16_var._rotation_is_identity) {
        if(!_Analyzer_16_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_16_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_16_var._position_relative, _Analyzer_16_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 21) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_16_var._position_relative, _Analyzer_16_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_16_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_16_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_16 [21] */
    /* begin component Analyzer_17=Monochromator_bent() [22] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_17_var._rotation_is_identity) {
        if(!_Analyzer_17_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_17_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_17_var._position_relative, _Analyzer_17_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 22) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_17_var._position_relative, _Analyzer_17_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_17_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_17_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_17 [22] */
    /* begin component Analyzer_18=Monochromator_bent() [23] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_18_var._rotation_is_identity) {
        if(!_Analyzer_18_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_18_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_18_var._position_relative, _Analyzer_18_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 23) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_18_var._position_relative, _Analyzer_18_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_18_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_18_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_18 [23] */
    /* begin component Analyzer_19=Monochromator_bent() [24] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_19_var._rotation_is_identity) {
        if(!_Analyzer_19_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_19_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_19_var._position_relative, _Analyzer_19_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 24) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_19_var._position_relative, _Analyzer_19_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_19_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_19_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_19 [24] */
    /* begin component Analyzer_20=Monochromator_bent() [25] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_20_var._rotation_is_identity) {
        if(!_Analyzer_20_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_20_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_20_var._position_relative, _Analyzer_20_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 25) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_20_var._position_relative, _Analyzer_20_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_20_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_20_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_20 [25] */
    /* begin component Analyzer_21=Monochromator_bent() [26] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_21_var._rotation_is_identity) {
        if(!_Analyzer_21_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_21_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_21_var._position_relative, _Analyzer_21_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 26) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_21_var._position_relative, _Analyzer_21_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_21_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_21_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_21 [26] */
    /* begin component Analyzer_22=Monochromator_bent() [27] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_22_var._rotation_is_identity) {
        if(!_Analyzer_22_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_22_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_22_var._position_relative, _Analyzer_22_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 27) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_22_var._position_relative, _Analyzer_22_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_22_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_22_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_22 [27] */
    /* begin component Analyzer_23=Monochromator_bent() [28] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_23_var._rotation_is_identity) {
        if(!_Analyzer_23_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_23_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_23_var._position_relative, _Analyzer_23_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 28) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_23_var._position_relative, _Analyzer_23_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_23_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_23_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_23 [28] */
    /* begin component Analyzer_24=Monochromator_bent() [29] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_24_var._rotation_is_identity) {
        if(!_Analyzer_24_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_24_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_24_var._position_relative, _Analyzer_24_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 29) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_24_var._position_relative, _Analyzer_24_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_24_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_24_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_24 [29] */
    /* begin component Analyzer_25=Monochromator_bent() [30] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_25_var._rotation_is_identity) {
        if(!_Analyzer_25_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_25_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_25_var._position_relative, _Analyzer_25_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 30) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_25_var._position_relative, _Analyzer_25_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_25_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_25_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_25 [30] */
    /* begin component Analyzer_26=Monochromator_bent() [31] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_26_var._rotation_is_identity) {
        if(!_Analyzer_26_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_26_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_26_var._position_relative, _Analyzer_26_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 31) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_26_var._position_relative, _Analyzer_26_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_26_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_26_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_26 [31] */
    /* begin component Analyzer_27=Monochromator_bent() [32] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_27_var._rotation_is_identity) {
        if(!_Analyzer_27_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_27_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_27_var._position_relative, _Analyzer_27_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 32) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_27_var._position_relative, _Analyzer_27_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_27_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_27_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_27 [32] */
    /* begin component Analyzer_28=Monochromator_bent() [33] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_28_var._rotation_is_identity) {
        if(!_Analyzer_28_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_28_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_28_var._position_relative, _Analyzer_28_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 33) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_28_var._position_relative, _Analyzer_28_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_28_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_28_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_28 [33] */
    /* begin component Analyzer_29=Monochromator_bent() [34] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_29_var._rotation_is_identity) {
        if(!_Analyzer_29_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_29_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_29_var._position_relative, _Analyzer_29_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 34) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_29_var._position_relative, _Analyzer_29_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_29_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_29_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_29 [34] */
    /* begin component Analyzer_30=Monochromator_bent() [35] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_30_var._rotation_is_identity) {
        if(!_Analyzer_30_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_30_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_30_var._position_relative, _Analyzer_30_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 35) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_30_var._position_relative, _Analyzer_30_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_30_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_30_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_30 [35] */
    /* begin component Analyzer_31=Monochromator_bent() [36] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_31_var._rotation_is_identity) {
        if(!_Analyzer_31_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_31_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_31_var._position_relative, _Analyzer_31_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 36) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_31_var._position_relative, _Analyzer_31_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_31_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_31_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_31 [36] */
    /* begin component Analyzer_32=Monochromator_bent() [37] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_32_var._rotation_is_identity) {
        if(!_Analyzer_32_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_32_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_32_var._position_relative, _Analyzer_32_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 37) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_32_var._position_relative, _Analyzer_32_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_32_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_32_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_32 [37] */
    /* begin component Analyzer_33=Monochromator_bent() [38] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_33_var._rotation_is_identity) {
        if(!_Analyzer_33_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_33_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_33_var._position_relative, _Analyzer_33_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 38) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_33_var._position_relative, _Analyzer_33_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_33_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_33_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_33 [38] */
    /* begin component Analyzer_34=Monochromator_bent() [39] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_34_var._rotation_is_identity) {
        if(!_Analyzer_34_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_34_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_34_var._position_relative, _Analyzer_34_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 39) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_34_var._position_relative, _Analyzer_34_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_34_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_34_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_34 [39] */
    /* begin component Analyzer_35=Monochromator_bent() [40] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_35_var._rotation_is_identity) {
        if(!_Analyzer_35_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_35_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_35_var._position_relative, _Analyzer_35_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 40) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_35_var._position_relative, _Analyzer_35_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_35_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_35_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_35 [40] */
    /* begin component Analyzer_36=Monochromator_bent() [41] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_36_var._rotation_is_identity) {
        if(!_Analyzer_36_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_36_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_36_var._position_relative, _Analyzer_36_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 41) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_36_var._position_relative, _Analyzer_36_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_36_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_36_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_36 [41] */
    /* begin component Analyzer_37=Monochromator_bent() [42] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_37_var._rotation_is_identity) {
        if(!_Analyzer_37_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_37_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_37_var._position_relative, _Analyzer_37_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 42) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_37_var._position_relative, _Analyzer_37_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_37_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_37_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_37 [42] */
    /* begin component Analyzer_38=Monochromator_bent() [43] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_38_var._rotation_is_identity) {
        if(!_Analyzer_38_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_38_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_38_var._position_relative, _Analyzer_38_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 43) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_38_var._position_relative, _Analyzer_38_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_38_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_38_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_38 [43] */
    /* begin component Analyzer_39=Monochromator_bent() [44] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_39_var._rotation_is_identity) {
        if(!_Analyzer_39_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_39_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_39_var._position_relative, _Analyzer_39_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 44) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_39_var._position_relative, _Analyzer_39_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_39_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_39_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_39 [44] */
    /* begin component Analyzer_40=Monochromator_bent() [45] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_40_var._rotation_is_identity) {
        if(!_Analyzer_40_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_40_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_40_var._position_relative, _Analyzer_40_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 45) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_40_var._position_relative, _Analyzer_40_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_40_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_40_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_40 [45] */
    /* begin component Analyzer_41=Monochromator_bent() [46] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_41_var._rotation_is_identity) {
        if(!_Analyzer_41_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_41_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_41_var._position_relative, _Analyzer_41_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 46) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_41_var._position_relative, _Analyzer_41_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_41_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_41_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_41 [46] */
    /* begin component Analyzer_42=Monochromator_bent() [47] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_42_var._rotation_is_identity) {
        if(!_Analyzer_42_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_42_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_42_var._position_relative, _Analyzer_42_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 47) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_42_var._position_relative, _Analyzer_42_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_42_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_42_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_42 [47] */
    /* begin component Analyzer_43=Monochromator_bent() [48] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_43_var._rotation_is_identity) {
        if(!_Analyzer_43_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_43_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_43_var._position_relative, _Analyzer_43_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 48) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_43_var._position_relative, _Analyzer_43_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_43_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_43_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_43 [48] */
    /* begin component Analyzer_44=Monochromator_bent() [49] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_44_var._rotation_is_identity) {
        if(!_Analyzer_44_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_44_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_44_var._position_relative, _Analyzer_44_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 49) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_44_var._position_relative, _Analyzer_44_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_44_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_44_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else particle_restore(_particle, &_particle_save); // not SCATTERED in GROUP, restore
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_44 [49] */
    /* begin component Analyzer_45=Monochromator_bent() [50] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Analyzer_45_var._rotation_is_identity) {
        if(!_Analyzer_45_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_45_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Analyzer_45_var._position_relative, _Analyzer_45_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 50) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      // 2nd or higher GROUP member, "reuse" coordinate-changed _particle_save from 1st GROUP element.
      mccoordschange(_Analyzer_45_var._position_relative, _Analyzer_45_var._rotation_relative, &_particle_save);
      DEBUG_COMP(_Analyzer_45_var._name);
      DEBUG_STATE();
      class_Monochromator_bent_trace(&_Analyzer_45_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORB;     // not SCATTERED at end of GROUP: removes left events
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Analyzer_45 [50] */
    /* begin component rotator=Arm() [51] */
    if (!ABSORBED && _particle->_index == 51) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component rotator [51] */
    /* begin component detector_pos=Arm() [52] */
    if (!ABSORBED && _particle->_index == 52) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component detector_pos [52] */
    /* begin component detector=Arm() [53] */
    if (!ABSORBED && _particle->_index == 53) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component detector [53] */
    /* begin component psd_monitor_end=PSD_monitor() [54] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_psd_monitor_end_var._rotation_is_identity) {
        if(!_psd_monitor_end_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _psd_monitor_end_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_psd_monitor_end_var._position_relative, _psd_monitor_end_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 54) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_psd_monitor_end_var._name);
      DEBUG_STATE();
      if ((( scat == 1 ))) // conditional WHEN execution
      class_PSD_monitor_trace(&_psd_monitor_end_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component psd_monitor_end [54] */
    /* begin component E_PSD_mon_end=Monitor_nD() [55] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_E_PSD_mon_end_var._rotation_is_identity) {
        if(!_E_PSD_mon_end_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _E_PSD_mon_end_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_E_PSD_mon_end_var._position_relative, _E_PSD_mon_end_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 55) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_E_PSD_mon_end_var._name);
      DEBUG_STATE();
      if ((( scat == 1 ))) // conditional WHEN execution
      class_Monitor_nD_trace(&_E_PSD_mon_end_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component E_PSD_mon_end [55] */
    if (_particle->_index > 55)
      ABSORBED++; /* absorbed when passed all components */
  } /* while !ABSORBED */

  DEBUG_LEAVE()
  particle_restore(_particle, &_particle_save);
  DEBUG_STATE()

  return(_particle->_index);
} /* raytrace */

/* loop to generate events and call raytrace() propagate them */
void raytrace_all(unsigned long long ncount, unsigned long seed) {

  /* CPU-loop */
  unsigned long long loops;
  loops = ceil((double)ncount/gpu_innerloop);
  /* if on GPU, printf has been globally nullified, re-enable here */
  #ifdef OPENACC
  #ifndef MULTICORE
  #undef printf
  #endif
  #endif

  #ifdef OPENACC
  if (ncount>gpu_innerloop) {
    printf("Defining %llu CPU loops around GPU kernel and adjusting ncount\n",loops);
    mcset_ncount(loops*gpu_innerloop);
  } else {
    #endif
    loops=1;
    gpu_innerloop = ncount;
    #ifdef OPENACC
  }
    #endif

  for (unsigned long long cloop=0; cloop<loops; cloop++) {
    #ifdef OPENACC
    if (loops>1) fprintf(stdout, "%d..", (int)cloop); fflush(stdout);
    #endif

    /* if on GPU, re-nullify printf */
    #ifdef OPENACC
    #ifndef MULTICORE
    #define printf(...) noprintf()
    #endif
    #endif

    #pragma acc parallel loop num_gangs(numgangs) vector_length(vecsize)
    for (unsigned long pidx=0 ; pidx < gpu_innerloop ; pidx++) {
      _class_particle particleN = mcgenstate(); // initial particle
      _class_particle* _particle = &particleN;
      particleN._uid = pidx;

      srandom(_hash((pidx+1)*(seed+1)));
      particle_uservar_init(_particle);

      raytrace(_particle);
    } /* inner for */
    seed = seed+gpu_innerloop;
  } /* CPU for */
  /* if on GPU, printf has been globally nullified, re-enable here */
  #ifdef OPENACC
  #ifndef MULTICORE
  #undef printf
  #endif
  #endif
  MPI_MASTER(
  printf("*** TRACE end *** \n");
  );
} /* raytrace_all */

#endif //no-FUNNEL

#ifdef FUNNEL
// Alternative raytrace algorithm which iterates all particles through
// one component at the time, can remove absorbs from the next loop and
// switch between cpu/gpu.
void raytrace_all_funnel(unsigned long long ncount, unsigned long seed) {

  // set up outer (CPU) loop / particle batches
  unsigned long long loops;

  /* if on GPU, printf has been globally nullified, re-enable here */
  #ifdef OPENACC
  #ifndef MULTICORE
  #undef printf
  #endif
  #endif

  #ifdef OPENACC
  loops = ceil((double)ncount/gpu_innerloop);
  if (ncount>gpu_innerloop) {
    printf("Defining %llu CPU loops around kernel and adjusting ncount\n",loops);
    mcset_ncount(loops*gpu_innerloop);
  } else {
  #endif
    loops=1;
    gpu_innerloop = ncount;
  #ifdef OPENACC
  }
  #endif

  // create particles struct and pointer arrays (same memory used by all batches)
  _class_particle* particles = malloc(gpu_innerloop*sizeof(_class_particle));
  _class_particle* pbuffer = malloc(gpu_innerloop*sizeof(_class_particle));
  long livebatchsize = gpu_innerloop;

  #undef ABSORB0
  #undef ABSORB
  #define ABSORB0 do { DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; } while(0)
  #define ABSORB ABSORB0
  // outer loop / particle batches
  for (unsigned long long cloop=0; cloop<loops; cloop++) {
    if (loops>1) fprintf(stdout, "%d..", (int)cloop); fflush(stdout);

    // init particles
    #pragma acc parallel loop present(particles)
    for (unsigned long pidx=0 ; pidx < livebatchsize ; pidx++) {
      // generate particle state, set loop index and seed
      particles[pidx] = mcgenstate();
      _class_particle* _particle = particles + pidx;
      _particle->_uid = pidx;
      srandom(_hash((pidx+1)*(seed+1))); // _particle->state usage built into srandom macro
      particle_uservar_init(_particle);
    }

    // iterate components

    #pragma acc parallel loop present(particles)
    for (unsigned long pidx=0 ; pidx < livebatchsize ; pidx++) {
      _class_particle* _particle = &particles[pidx];
      _class_particle _particle_save;

      // origin
    if (!ABSORBED && _particle->_index == 1) {
        if (_origin_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _origin_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_origin_var._position_relative, _origin_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Progress_bar_trace(&_origin_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // source
    if (!ABSORBED && _particle->_index == 2) {
        if (_source_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _source_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_source_var._position_relative, _source_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Source_div_trace(&_source_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // psd_test
    if (!ABSORBED && _particle->_index == 3) {
        if (_psd_test_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _psd_test_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_psd_test_var._position_relative, _psd_test_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_PSD_monitor_trace(&_psd_test_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // e_monitor
    if (!ABSORBED && _particle->_index == 4) {
        if (_e_monitor_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _e_monitor_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_e_monitor_var._position_relative, _e_monitor_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_E_monitor_trace(&_e_monitor_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // Analyzer_0
    if (!ABSORBED && _particle->_index == 5) {
        if (_Analyzer_0_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_0_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_0_var._position_relative, _Analyzer_0_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_0_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_1
    if (!ABSORBED && _particle->_index == 6) {
        if (_Analyzer_1_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_1_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_1_var._position_relative, _Analyzer_1_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_1_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_2
    if (!ABSORBED && _particle->_index == 7) {
        if (_Analyzer_2_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_2_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_2_var._position_relative, _Analyzer_2_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_2_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_3
    if (!ABSORBED && _particle->_index == 8) {
        if (_Analyzer_3_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_3_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_3_var._position_relative, _Analyzer_3_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_3_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_4
    if (!ABSORBED && _particle->_index == 9) {
        if (_Analyzer_4_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_4_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_4_var._position_relative, _Analyzer_4_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_4_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_5
    if (!ABSORBED && _particle->_index == 10) {
        if (_Analyzer_5_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_5_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_5_var._position_relative, _Analyzer_5_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_5_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_6
    if (!ABSORBED && _particle->_index == 11) {
        if (_Analyzer_6_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_6_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_6_var._position_relative, _Analyzer_6_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_6_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_7
    if (!ABSORBED && _particle->_index == 12) {
        if (_Analyzer_7_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_7_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_7_var._position_relative, _Analyzer_7_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_7_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_8
    if (!ABSORBED && _particle->_index == 13) {
        if (_Analyzer_8_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_8_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_8_var._position_relative, _Analyzer_8_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_8_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_9
    if (!ABSORBED && _particle->_index == 14) {
        if (_Analyzer_9_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_9_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_9_var._position_relative, _Analyzer_9_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_9_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_10
    if (!ABSORBED && _particle->_index == 15) {
        if (_Analyzer_10_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_10_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_10_var._position_relative, _Analyzer_10_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_10_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_11
    if (!ABSORBED && _particle->_index == 16) {
        if (_Analyzer_11_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_11_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_11_var._position_relative, _Analyzer_11_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_11_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_12
    if (!ABSORBED && _particle->_index == 17) {
        if (_Analyzer_12_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_12_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_12_var._position_relative, _Analyzer_12_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_12_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_13
    if (!ABSORBED && _particle->_index == 18) {
        if (_Analyzer_13_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_13_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_13_var._position_relative, _Analyzer_13_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_13_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_14
    if (!ABSORBED && _particle->_index == 19) {
        if (_Analyzer_14_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_14_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_14_var._position_relative, _Analyzer_14_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_14_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_15
    if (!ABSORBED && _particle->_index == 20) {
        if (_Analyzer_15_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_15_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_15_var._position_relative, _Analyzer_15_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_15_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_16
    if (!ABSORBED && _particle->_index == 21) {
        if (_Analyzer_16_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_16_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_16_var._position_relative, _Analyzer_16_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_16_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_17
    if (!ABSORBED && _particle->_index == 22) {
        if (_Analyzer_17_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_17_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_17_var._position_relative, _Analyzer_17_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_17_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_18
    if (!ABSORBED && _particle->_index == 23) {
        if (_Analyzer_18_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_18_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_18_var._position_relative, _Analyzer_18_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_18_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_19
    if (!ABSORBED && _particle->_index == 24) {
        if (_Analyzer_19_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_19_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_19_var._position_relative, _Analyzer_19_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_19_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_20
    if (!ABSORBED && _particle->_index == 25) {
        if (_Analyzer_20_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_20_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_20_var._position_relative, _Analyzer_20_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_20_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_21
    if (!ABSORBED && _particle->_index == 26) {
        if (_Analyzer_21_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_21_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_21_var._position_relative, _Analyzer_21_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_21_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_22
    if (!ABSORBED && _particle->_index == 27) {
        if (_Analyzer_22_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_22_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_22_var._position_relative, _Analyzer_22_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_22_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_23
    if (!ABSORBED && _particle->_index == 28) {
        if (_Analyzer_23_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_23_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_23_var._position_relative, _Analyzer_23_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_23_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_24
    if (!ABSORBED && _particle->_index == 29) {
        if (_Analyzer_24_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_24_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_24_var._position_relative, _Analyzer_24_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_24_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_25
    if (!ABSORBED && _particle->_index == 30) {
        if (_Analyzer_25_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_25_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_25_var._position_relative, _Analyzer_25_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_25_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_26
    if (!ABSORBED && _particle->_index == 31) {
        if (_Analyzer_26_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_26_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_26_var._position_relative, _Analyzer_26_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_26_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_27
    if (!ABSORBED && _particle->_index == 32) {
        if (_Analyzer_27_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_27_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_27_var._position_relative, _Analyzer_27_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_27_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_28
    if (!ABSORBED && _particle->_index == 33) {
        if (_Analyzer_28_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_28_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_28_var._position_relative, _Analyzer_28_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_28_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_29
    if (!ABSORBED && _particle->_index == 34) {
        if (_Analyzer_29_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_29_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_29_var._position_relative, _Analyzer_29_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_29_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_30
    if (!ABSORBED && _particle->_index == 35) {
        if (_Analyzer_30_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_30_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_30_var._position_relative, _Analyzer_30_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_30_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_31
    if (!ABSORBED && _particle->_index == 36) {
        if (_Analyzer_31_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_31_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_31_var._position_relative, _Analyzer_31_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_31_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_32
    if (!ABSORBED && _particle->_index == 37) {
        if (_Analyzer_32_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_32_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_32_var._position_relative, _Analyzer_32_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_32_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_33
    if (!ABSORBED && _particle->_index == 38) {
        if (_Analyzer_33_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_33_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_33_var._position_relative, _Analyzer_33_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_33_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_34
    if (!ABSORBED && _particle->_index == 39) {
        if (_Analyzer_34_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_34_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_34_var._position_relative, _Analyzer_34_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_34_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_35
    if (!ABSORBED && _particle->_index == 40) {
        if (_Analyzer_35_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_35_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_35_var._position_relative, _Analyzer_35_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_35_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_36
    if (!ABSORBED && _particle->_index == 41) {
        if (_Analyzer_36_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_36_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_36_var._position_relative, _Analyzer_36_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_36_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_37
    if (!ABSORBED && _particle->_index == 42) {
        if (_Analyzer_37_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_37_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_37_var._position_relative, _Analyzer_37_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_37_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_38
    if (!ABSORBED && _particle->_index == 43) {
        if (_Analyzer_38_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_38_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_38_var._position_relative, _Analyzer_38_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_38_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_39
    if (!ABSORBED && _particle->_index == 44) {
        if (_Analyzer_39_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_39_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_39_var._position_relative, _Analyzer_39_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_39_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_40
    if (!ABSORBED && _particle->_index == 45) {
        if (_Analyzer_40_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_40_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_40_var._position_relative, _Analyzer_40_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_40_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_41
    if (!ABSORBED && _particle->_index == 46) {
        if (_Analyzer_41_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_41_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_41_var._position_relative, _Analyzer_41_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_41_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_42
    if (!ABSORBED && _particle->_index == 47) {
        if (_Analyzer_42_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_42_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_42_var._position_relative, _Analyzer_42_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_42_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_43
    if (!ABSORBED && _particle->_index == 48) {
        if (_Analyzer_43_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_43_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_43_var._position_relative, _Analyzer_43_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_43_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_44
    if (!ABSORBED && _particle->_index == 49) {
        if (_Analyzer_44_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_44_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_44_var._position_relative, _Analyzer_44_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_44_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORBED=0; // not SCATTERED within GROUP: always tries next
        _particle->_index++;
      }

      // Analyzer_45
    if (!ABSORBED && _particle->_index == 50) {
        if (_Analyzer_45_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Analyzer_45_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Analyzer_45_var._position_relative, _Analyzer_45_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_bent_trace(&_Analyzer_45_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      // GROUP Analyzer: from Analyzer_0 [5] to Analyzer_45 [50]
      if (SCATTERED) _particle->_index = 50; // when SCATTERED in GROUP: reach exit of GROUP after Analyzer_45
      else ABSORB;     // not SCATTERED at end of GROUP: removes left events
        _particle->_index++;
      }

      // rotator
    if (!ABSORBED && _particle->_index == 51) {
        _particle->_index++;
      }

      // detector_pos
    if (!ABSORBED && _particle->_index == 52) {
        _particle->_index++;
      }

      // detector
    if (!ABSORBED && _particle->_index == 53) {
        _particle->_index++;
      }

      // psd_monitor_end
    if (!ABSORBED && _particle->_index == 54) {
        if (_psd_monitor_end_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _psd_monitor_end_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_psd_monitor_end_var._position_relative, _psd_monitor_end_var._rotation_relative, _particle);
        _particle_save = *_particle;
        if ((( scat == 1 ))) // conditional WHEN
        class_PSD_monitor_trace(&_psd_monitor_end_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // E_PSD_mon_end
    if (!ABSORBED && _particle->_index == 55) {
        if (_E_PSD_mon_end_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _E_PSD_mon_end_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_E_PSD_mon_end_var._position_relative, _E_PSD_mon_end_var._rotation_relative, _particle);
        _particle_save = *_particle;
        if ((( scat == 1 ))) // conditional WHEN
        class_Monitor_nD_trace(&_E_PSD_mon_end_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

    }

    // jump to next viable seed
    seed = seed + gpu_innerloop;
  } // outer loop / particle batches

  free(particles);
  free(pbuffer);

  printf("\n");
} /* raytrace_all_funnel */
#endif // FUNNEL

#undef x
#undef y
#undef z
#undef vx
#undef vy
#undef vz
#undef t
#undef sx
#undef sy
#undef sz
#undef p
#undef mcgravitation
#undef mcMagnet
#undef allow_backprop
#undef _mctmp_a
#undef _mctmp_b
#undef _mctmp_c
#ifdef OPENACC
#ifndef MULTICORE
#undef strlen
#undef strcmp
#undef exit
#undef printf
#undef sprintf
#undef fprintf
#endif
#endif
#undef SCATTERED
#undef RESTORE
#undef RESTORE_NEUTRON
#undef STORE_NEUTRON
#undef ABSORBED
#undef ABSORB
#undef ABSORB0
/* *****************************************************************************
* instrument 'template_simple' and components SAVE
***************************************************************************** */

_class_Progress_bar *class_Progress_bar_save(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_save] component origin=Progress_bar() SAVE [C:\\mcstas-3.4\\lib\\misc\\Progress_bar.comp:120]");

  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", instrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, instrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &(instrument->counter_N[1]),&(instrument->counter_P[1]),&(instrument->counter_P2[1]),
        filename);

  }
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_save */

_class_PSD_monitor *class_PSD_monitor_save(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_psd_test_save] component psd_test=PSD_monitor() SAVE [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:106]");

    if (!nowritefile) {
      DETECTOR_OUT_2D(
          "PSD monitor",
          "X position [cm]",
          "Y position [cm]",
          xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
          nx, ny,
          &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
          filename);
    }
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_save */

_class_E_monitor *class_E_monitor_save(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_save] component e_monitor=E_monitor() SAVE [C:\\mcstas-3.4\\lib\\monitors\\E_monitor.comp:131]");

if (!nowritefile) {
  DETECTOR_OUT_1D(
      "Energy monitor",
      "Energy [meV]",
      "Intensity",
      "E", Emin, Emax, nE,
      &E_N[0],&E_p[0],&E_p2[0],
      filename);
  if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
   S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
}
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_save */

_class_Monitor_nD *class_Monitor_nD_save(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  SIG_MESSAGE("[_E_PSD_mon_end_save] component E_PSD_mon_end=Monitor_nD() SAVE [C:\\mcstas-3.4\\lib\\monitors\\Monitor_nD.comp:615]");

if (!nowritefile) {
  /* save results, but do not free pointers */
  detector = Monitor_nD_Save(&DEFS, &Vars);
}
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef nowritefile
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_save */



int save(FILE *handle) { /* called by mccode_main for template_simple:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(&_origin_var);


  class_PSD_monitor_save(&_psd_test_var);

  class_E_monitor_save(&_e_monitor_var);


















































  class_PSD_monitor_save(&_psd_monitor_end_var);

  class_Monitor_nD_save(&_E_PSD_mon_end_var);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'template_simple' and components FINALLY
***************************************************************************** */

_class_Progress_bar *class_Progress_bar_finally(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_finally] component origin=Progress_bar() FINALLY [C:\\mcstas-3.4\\lib\\misc\\Progress_bar.comp:138]");

  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", instrument_name, dirname ? dirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3600.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_finally */

_class_PSD_monitor *class_PSD_monitor_finally(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_psd_test_finally] component psd_test=PSD_monitor() FINALLY [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:119]");

  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_finally */

_class_E_monitor *class_E_monitor_finally(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_finally] component e_monitor=E_monitor() FINALLY [C:\\mcstas-3.4\\lib\\monitors\\E_monitor.comp:146]");

  destroy_darr1d(E_N);
  destroy_darr1d(E_p);
  destroy_darr1d(E_p2);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_finally */

_class_Monochromator_bent *class_Monochromator_bent_finally(_class_Monochromator_bent *_comp
) {
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define xthickness (_comp->_parameters.xthickness)
  #define radius_x (_comp->_parameters.radius_x)
  #define plane_of_reflection (_comp->_parameters.plane_of_reflection)
  #define angle_to_cut_horizontal (_comp->_parameters.angle_to_cut_horizontal)
  #define angle_to_cut_vertical (_comp->_parameters.angle_to_cut_vertical)
  #define mosaicity (_comp->_parameters.mosaicity)
  #define domainthickness (_comp->_parameters.domainthickness)
  #define temperature (_comp->_parameters.temperature)
  #define verbose (_comp->_parameters.verbose)
  #define neutron (_comp->_parameters.neutron)
  #define monochromator (_comp->_parameters.monochromator)
  #define non_scattered (_comp->_parameters.non_scattered)
  #define scattered (_comp->_parameters.scattered)
  #define non_hit (_comp->_parameters.non_hit)
  SIG_MESSAGE("[_Analyzer_0_finally] component Analyzer_0=Monochromator_bent() FINALLY [Monochromator_bent.comp:957]");

	double beta = neutron.beta;
	double eps_zero = neutron.eps_zero;
	int scatter_count = scattered;
  #undef zwidth
  #undef yheight
  #undef xthickness
  #undef radius_x
  #undef plane_of_reflection
  #undef angle_to_cut_horizontal
  #undef angle_to_cut_vertical
  #undef mosaicity
  #undef domainthickness
  #undef temperature
  #undef verbose
  #undef neutron
  #undef monochromator
  #undef non_scattered
  #undef scattered
  #undef non_hit
  return(_comp);
} /* class_Monochromator_bent_finally */

_class_Monitor_nD *class_Monitor_nD_finally(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  SIG_MESSAGE("[_E_PSD_mon_end_finally] component E_PSD_mon_end=Monitor_nD() FINALLY [C:\\mcstas-3.4\\lib\\monitors\\Monitor_nD.comp:623]");

  /* free pointers */
  Monitor_nD_Finally(&DEFS, &Vars);
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef nowritefile
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_finally */



int finally(void) { /* called by mccode_main for template_simple:FINALLY */
#pragma acc update host(_origin_var)
#pragma acc update host(_source_var)
#pragma acc update host(_psd_test_var)
#pragma acc update host(_e_monitor_var)
#pragma acc update host(_Analyzer_0_var)
#pragma acc update host(_Analyzer_1_var)
#pragma acc update host(_Analyzer_2_var)
#pragma acc update host(_Analyzer_3_var)
#pragma acc update host(_Analyzer_4_var)
#pragma acc update host(_Analyzer_5_var)
#pragma acc update host(_Analyzer_6_var)
#pragma acc update host(_Analyzer_7_var)
#pragma acc update host(_Analyzer_8_var)
#pragma acc update host(_Analyzer_9_var)
#pragma acc update host(_Analyzer_10_var)
#pragma acc update host(_Analyzer_11_var)
#pragma acc update host(_Analyzer_12_var)
#pragma acc update host(_Analyzer_13_var)
#pragma acc update host(_Analyzer_14_var)
#pragma acc update host(_Analyzer_15_var)
#pragma acc update host(_Analyzer_16_var)
#pragma acc update host(_Analyzer_17_var)
#pragma acc update host(_Analyzer_18_var)
#pragma acc update host(_Analyzer_19_var)
#pragma acc update host(_Analyzer_20_var)
#pragma acc update host(_Analyzer_21_var)
#pragma acc update host(_Analyzer_22_var)
#pragma acc update host(_Analyzer_23_var)
#pragma acc update host(_Analyzer_24_var)
#pragma acc update host(_Analyzer_25_var)
#pragma acc update host(_Analyzer_26_var)
#pragma acc update host(_Analyzer_27_var)
#pragma acc update host(_Analyzer_28_var)
#pragma acc update host(_Analyzer_29_var)
#pragma acc update host(_Analyzer_30_var)
#pragma acc update host(_Analyzer_31_var)
#pragma acc update host(_Analyzer_32_var)
#pragma acc update host(_Analyzer_33_var)
#pragma acc update host(_Analyzer_34_var)
#pragma acc update host(_Analyzer_35_var)
#pragma acc update host(_Analyzer_36_var)
#pragma acc update host(_Analyzer_37_var)
#pragma acc update host(_Analyzer_38_var)
#pragma acc update host(_Analyzer_39_var)
#pragma acc update host(_Analyzer_40_var)
#pragma acc update host(_Analyzer_41_var)
#pragma acc update host(_Analyzer_42_var)
#pragma acc update host(_Analyzer_43_var)
#pragma acc update host(_Analyzer_44_var)
#pragma acc update host(_Analyzer_45_var)
#pragma acc update host(_rotator_var)
#pragma acc update host(_detector_pos_var)
#pragma acc update host(_detector_var)
#pragma acc update host(_psd_monitor_end_var)
#pragma acc update host(_E_PSD_mon_end_var)
#pragma acc update host(_instrument_var)

  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(&_origin_var);


  class_PSD_monitor_finally(&_psd_test_var);

  class_E_monitor_finally(&_e_monitor_var);

  class_Monochromator_bent_finally(&_Analyzer_0_var);

  class_Monochromator_bent_finally(&_Analyzer_1_var);

  class_Monochromator_bent_finally(&_Analyzer_2_var);

  class_Monochromator_bent_finally(&_Analyzer_3_var);

  class_Monochromator_bent_finally(&_Analyzer_4_var);

  class_Monochromator_bent_finally(&_Analyzer_5_var);

  class_Monochromator_bent_finally(&_Analyzer_6_var);

  class_Monochromator_bent_finally(&_Analyzer_7_var);

  class_Monochromator_bent_finally(&_Analyzer_8_var);

  class_Monochromator_bent_finally(&_Analyzer_9_var);

  class_Monochromator_bent_finally(&_Analyzer_10_var);

  class_Monochromator_bent_finally(&_Analyzer_11_var);

  class_Monochromator_bent_finally(&_Analyzer_12_var);

  class_Monochromator_bent_finally(&_Analyzer_13_var);

  class_Monochromator_bent_finally(&_Analyzer_14_var);

  class_Monochromator_bent_finally(&_Analyzer_15_var);

  class_Monochromator_bent_finally(&_Analyzer_16_var);

  class_Monochromator_bent_finally(&_Analyzer_17_var);

  class_Monochromator_bent_finally(&_Analyzer_18_var);

  class_Monochromator_bent_finally(&_Analyzer_19_var);

  class_Monochromator_bent_finally(&_Analyzer_20_var);

  class_Monochromator_bent_finally(&_Analyzer_21_var);

  class_Monochromator_bent_finally(&_Analyzer_22_var);

  class_Monochromator_bent_finally(&_Analyzer_23_var);

  class_Monochromator_bent_finally(&_Analyzer_24_var);

  class_Monochromator_bent_finally(&_Analyzer_25_var);

  class_Monochromator_bent_finally(&_Analyzer_26_var);

  class_Monochromator_bent_finally(&_Analyzer_27_var);

  class_Monochromator_bent_finally(&_Analyzer_28_var);

  class_Monochromator_bent_finally(&_Analyzer_29_var);

  class_Monochromator_bent_finally(&_Analyzer_30_var);

  class_Monochromator_bent_finally(&_Analyzer_31_var);

  class_Monochromator_bent_finally(&_Analyzer_32_var);

  class_Monochromator_bent_finally(&_Analyzer_33_var);

  class_Monochromator_bent_finally(&_Analyzer_34_var);

  class_Monochromator_bent_finally(&_Analyzer_35_var);

  class_Monochromator_bent_finally(&_Analyzer_36_var);

  class_Monochromator_bent_finally(&_Analyzer_37_var);

  class_Monochromator_bent_finally(&_Analyzer_38_var);

  class_Monochromator_bent_finally(&_Analyzer_39_var);

  class_Monochromator_bent_finally(&_Analyzer_40_var);

  class_Monochromator_bent_finally(&_Analyzer_41_var);

  class_Monochromator_bent_finally(&_Analyzer_42_var);

  class_Monochromator_bent_finally(&_Analyzer_43_var);

  class_Monochromator_bent_finally(&_Analyzer_44_var);

  class_Monochromator_bent_finally(&_Analyzer_45_var);




  class_PSD_monitor_finally(&_psd_monitor_end_var);

  class_Monitor_nD_finally(&_E_PSD_mon_end_var);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'template_simple' and components DISPLAY
***************************************************************************** */

  #define magnify     mcdis_magnify
  #define line        mcdis_line
  #define dashed_line mcdis_dashed_line
  #define multiline   mcdis_multiline
  #define rectangle   mcdis_rectangle
  #define box         mcdis_box
  #define circle      mcdis_circle
  #define cylinder    mcdis_cylinder
  #define sphere      mcdis_sphere
_class_Progress_bar *class_Progress_bar_display(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_display] component origin=Progress_bar() DISPLAY [C:\\mcstas-3.4\\lib\\misc\\Progress_bar.comp:152]");

  printf("MCDISPLAY: component %s\n", _comp->_name);

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_display */

_class_Source_div *class_Source_div_display(_class_Source_div *_comp
) {
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define E0 (_comp->_parameters.E0)
  #define dE (_comp->_parameters.dE)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define gauss (_comp->_parameters.gauss)
  #define flux (_comp->_parameters.flux)
  #define sigmah (_comp->_parameters.sigmah)
  #define sigmav (_comp->_parameters.sigmav)
  #define p_init (_comp->_parameters.p_init)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  SIG_MESSAGE("[_source_display] component source=Source_div() DISPLAY [C:\\mcstas-3.4\\lib\\sources\\Source_div.comp:170]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  multiline(5, -xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0, -yheight/2.0, 0.0);
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
  #undef xwidth
  #undef yheight
  #undef focus_aw
  #undef focus_ah
  #undef E0
  #undef dE
  #undef lambda0
  #undef dlambda
  #undef gauss
  #undef flux
  #undef sigmah
  #undef sigmav
  #undef p_init
  #undef dist
  #undef focus_xw
  #undef focus_yh
  return(_comp);
} /* class_Source_div_display */

_class_PSD_monitor *class_PSD_monitor_display(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_psd_test_display] component psd_test=PSD_monitor() DISPLAY [C:\\mcstas-3.4\\lib\\monitors\\PSD_monitor.comp:126]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_display */

_class_E_monitor *class_E_monitor_display(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_display] component e_monitor=E_monitor() DISPLAY [C:\\mcstas-3.4\\lib\\monitors\\E_monitor.comp:153]");

  printf("MCDISPLAY: component %s\n", _comp->_name);

  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_display */

_class_Monochromator_bent *class_Monochromator_bent_display(_class_Monochromator_bent *_comp
) {
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define xthickness (_comp->_parameters.xthickness)
  #define radius_x (_comp->_parameters.radius_x)
  #define plane_of_reflection (_comp->_parameters.plane_of_reflection)
  #define angle_to_cut_horizontal (_comp->_parameters.angle_to_cut_horizontal)
  #define angle_to_cut_vertical (_comp->_parameters.angle_to_cut_vertical)
  #define mosaicity (_comp->_parameters.mosaicity)
  #define domainthickness (_comp->_parameters.domainthickness)
  #define temperature (_comp->_parameters.temperature)
  #define verbose (_comp->_parameters.verbose)
  #define neutron (_comp->_parameters.neutron)
  #define monochromator (_comp->_parameters.monochromator)
  #define non_scattered (_comp->_parameters.non_scattered)
  #define scattered (_comp->_parameters.scattered)
  #define non_hit (_comp->_parameters.non_hit)
  SIG_MESSAGE("[_Analyzer_0_display] component Analyzer_0=Monochromator_bent() DISPLAY [Monochromator_bent.comp:964]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
	double x_edges[2];
	double z_edges[2];
	double y_edges[2];
	double angle_range = monochromator.max_angle - monochromator.min_angle;
	if (radius_x == 0){
		radius_x = 10000000;
	}

	x_edges[0] = radius_x - (radius_x-xthickness/2)*cos(angle_range/2);
	x_edges[1] = radius_x - (radius_x+xthickness/2)*cos(angle_range/2);

	z_edges[0] = (radius_x-xthickness/2)*sin(angle_range/2);
	z_edges[1] = (radius_x+xthickness/2)*sin(angle_range/2);

	y_edges[0] = -yheight/2;
	y_edges[1] = yheight/2;

	double temp_x[4], temp_z[4];
	temp_x[0] = x_edges[0];
	temp_x[2] = x_edges[1];
	temp_z[0] = - z_edges[0];
	temp_z[2] = - z_edges[1];

	int n = 20;
	double delta_theta = angle_range/(n-1);
	for (int i=1; i<n; i++){
		// Inner upper
		temp_x[1] = radius_x - (radius_x-xthickness/2)*cos(angle_range/2-delta_theta*i);
		temp_x[3] = radius_x - (radius_x+xthickness/2)*cos(angle_range/2-delta_theta*i);
		temp_z[1] = -(radius_x-xthickness/2)*sin(angle_range/2-delta_theta*i);
		temp_z[3] = -(radius_x+xthickness/2)*sin(angle_range/2-delta_theta*i);

		line(temp_x[0],y_edges[0],temp_z[0],
			temp_x[1],y_edges[0],temp_z[1]);
		line(temp_x[0],y_edges[1],temp_z[0],
			temp_x[1],y_edges[1],temp_z[1]);
		line(temp_x[2],y_edges[0],temp_z[2],
			temp_x[3],y_edges[0],temp_z[3]);
		line(temp_x[2],y_edges[1],temp_z[2],
			temp_x[3],y_edges[1],temp_z[3]);

		// New starting values
		temp_x[0] = temp_x[1];
		temp_x[2] = temp_x[3];
		temp_z[0] = temp_z[1];
		temp_z[2] = temp_z[3];

	}

	line(x_edges[0],y_edges[0],z_edges[0],
		x_edges[0],y_edges[1],z_edges[0]);
	line(x_edges[1],y_edges[0],z_edges[1],
		x_edges[1],y_edges[1],z_edges[1]);
	line(x_edges[0],y_edges[0],-z_edges[0],
		x_edges[0],y_edges[1],-z_edges[0]);
	line(x_edges[1],y_edges[0],-z_edges[1],
		x_edges[1],y_edges[1],-z_edges[1]);

	line(x_edges[0],y_edges[0],z_edges[0],
		x_edges[1],y_edges[0],z_edges[1]);
	line(x_edges[0],y_edges[1],z_edges[0],
		x_edges[1],y_edges[1],z_edges[1]);
	line(x_edges[0],y_edges[0],-z_edges[0],
		x_edges[1],y_edges[0],-z_edges[1]);
	line(x_edges[0],y_edges[1],-z_edges[0],
		x_edges[1],y_edges[1],-z_edges[1]);
  #undef zwidth
  #undef yheight
  #undef xthickness
  #undef radius_x
  #undef plane_of_reflection
  #undef angle_to_cut_horizontal
  #undef angle_to_cut_vertical
  #undef mosaicity
  #undef domainthickness
  #undef temperature
  #undef verbose
  #undef neutron
  #undef monochromator
  #undef non_scattered
  #undef scattered
  #undef non_hit
  return(_comp);
} /* class_Monochromator_bent_display */

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  SIG_MESSAGE("[_rotator_display] component rotator=Arm() DISPLAY [C:\\mcstas-3.4\\lib\\optics\\Arm.comp:40]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_Monitor_nD *class_Monitor_nD_display(_class_Monitor_nD *_comp
) {
  #define user1 (_comp->_parameters.user1)
  #define user2 (_comp->_parameters.user2)
  #define user3 (_comp->_parameters.user3)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define zdepth (_comp->_parameters.zdepth)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define bins (_comp->_parameters.bins)
  #define min (_comp->_parameters.min)
  #define max (_comp->_parameters.max)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define radius (_comp->_parameters.radius)
  #define options (_comp->_parameters.options)
  #define filename (_comp->_parameters.filename)
  #define geometry (_comp->_parameters.geometry)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define username1 (_comp->_parameters.username1)
  #define username2 (_comp->_parameters.username2)
  #define username3 (_comp->_parameters.username3)
  #define DEFS (_comp->_parameters.DEFS)
  #define Vars (_comp->_parameters.Vars)
  #define detector (_comp->_parameters.detector)
  #define offdata (_comp->_parameters.offdata)
  SIG_MESSAGE("[_E_PSD_mon_end_display] component E_PSD_mon_end=Monitor_nD() DISPLAY [C:\\mcstas-3.4\\lib\\monitors\\Monitor_nD.comp:629]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  if (geometry && strlen(geometry) && strcmp(geometry,"0") && strcmp(geometry, "NULL"))
  {
    off_display(offdata);
  } else {
    Monitor_nD_McDisplay(&DEFS, &Vars);
  }
  #undef user1
  #undef user2
  #undef user3
  #undef xwidth
  #undef yheight
  #undef zdepth
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef zmin
  #undef zmax
  #undef bins
  #undef min
  #undef max
  #undef restore_neutron
  #undef radius
  #undef options
  #undef filename
  #undef geometry
  #undef nowritefile
  #undef username1
  #undef username2
  #undef username3
  #undef DEFS
  #undef Vars
  #undef detector
  #undef offdata
  return(_comp);
} /* class_Monitor_nD_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for template_simple:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(&_origin_var);

  class_Source_div_display(&_source_var);

  class_PSD_monitor_display(&_psd_test_var);

  class_E_monitor_display(&_e_monitor_var);

  class_Monochromator_bent_display(&_Analyzer_0_var);

  class_Monochromator_bent_display(&_Analyzer_1_var);

  class_Monochromator_bent_display(&_Analyzer_2_var);

  class_Monochromator_bent_display(&_Analyzer_3_var);

  class_Monochromator_bent_display(&_Analyzer_4_var);

  class_Monochromator_bent_display(&_Analyzer_5_var);

  class_Monochromator_bent_display(&_Analyzer_6_var);

  class_Monochromator_bent_display(&_Analyzer_7_var);

  class_Monochromator_bent_display(&_Analyzer_8_var);

  class_Monochromator_bent_display(&_Analyzer_9_var);

  class_Monochromator_bent_display(&_Analyzer_10_var);

  class_Monochromator_bent_display(&_Analyzer_11_var);

  class_Monochromator_bent_display(&_Analyzer_12_var);

  class_Monochromator_bent_display(&_Analyzer_13_var);

  class_Monochromator_bent_display(&_Analyzer_14_var);

  class_Monochromator_bent_display(&_Analyzer_15_var);

  class_Monochromator_bent_display(&_Analyzer_16_var);

  class_Monochromator_bent_display(&_Analyzer_17_var);

  class_Monochromator_bent_display(&_Analyzer_18_var);

  class_Monochromator_bent_display(&_Analyzer_19_var);

  class_Monochromator_bent_display(&_Analyzer_20_var);

  class_Monochromator_bent_display(&_Analyzer_21_var);

  class_Monochromator_bent_display(&_Analyzer_22_var);

  class_Monochromator_bent_display(&_Analyzer_23_var);

  class_Monochromator_bent_display(&_Analyzer_24_var);

  class_Monochromator_bent_display(&_Analyzer_25_var);

  class_Monochromator_bent_display(&_Analyzer_26_var);

  class_Monochromator_bent_display(&_Analyzer_27_var);

  class_Monochromator_bent_display(&_Analyzer_28_var);

  class_Monochromator_bent_display(&_Analyzer_29_var);

  class_Monochromator_bent_display(&_Analyzer_30_var);

  class_Monochromator_bent_display(&_Analyzer_31_var);

  class_Monochromator_bent_display(&_Analyzer_32_var);

  class_Monochromator_bent_display(&_Analyzer_33_var);

  class_Monochromator_bent_display(&_Analyzer_34_var);

  class_Monochromator_bent_display(&_Analyzer_35_var);

  class_Monochromator_bent_display(&_Analyzer_36_var);

  class_Monochromator_bent_display(&_Analyzer_37_var);

  class_Monochromator_bent_display(&_Analyzer_38_var);

  class_Monochromator_bent_display(&_Analyzer_39_var);

  class_Monochromator_bent_display(&_Analyzer_40_var);

  class_Monochromator_bent_display(&_Analyzer_41_var);

  class_Monochromator_bent_display(&_Analyzer_42_var);

  class_Monochromator_bent_display(&_Analyzer_43_var);

  class_Monochromator_bent_display(&_Analyzer_44_var);

  class_Monochromator_bent_display(&_Analyzer_45_var);

  class_Arm_display(&_rotator_var);

  class_Arm_display(&_detector_pos_var);

  class_Arm_display(&_detector_var);

  class_PSD_monitor_display(&_psd_monitor_end_var);

  class_Monitor_nD_display(&_E_PSD_mon_end_var);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  #ifdef OPENACC
    #define strcmp(a,b) str_comp(a,b)
  #endif
  if (!strcmp(compname, "origin")) return (void *) &(_origin_var._parameters);
  if (!strcmp(compname, "source")) return (void *) &(_source_var._parameters);
  if (!strcmp(compname, "psd_test")) return (void *) &(_psd_test_var._parameters);
  if (!strcmp(compname, "e_monitor")) return (void *) &(_e_monitor_var._parameters);
  if (!strcmp(compname, "Analyzer_0")) return (void *) &(_Analyzer_0_var._parameters);
  if (!strcmp(compname, "Analyzer_1")) return (void *) &(_Analyzer_1_var._parameters);
  if (!strcmp(compname, "Analyzer_2")) return (void *) &(_Analyzer_2_var._parameters);
  if (!strcmp(compname, "Analyzer_3")) return (void *) &(_Analyzer_3_var._parameters);
  if (!strcmp(compname, "Analyzer_4")) return (void *) &(_Analyzer_4_var._parameters);
  if (!strcmp(compname, "Analyzer_5")) return (void *) &(_Analyzer_5_var._parameters);
  if (!strcmp(compname, "Analyzer_6")) return (void *) &(_Analyzer_6_var._parameters);
  if (!strcmp(compname, "Analyzer_7")) return (void *) &(_Analyzer_7_var._parameters);
  if (!strcmp(compname, "Analyzer_8")) return (void *) &(_Analyzer_8_var._parameters);
  if (!strcmp(compname, "Analyzer_9")) return (void *) &(_Analyzer_9_var._parameters);
  if (!strcmp(compname, "Analyzer_10")) return (void *) &(_Analyzer_10_var._parameters);
  if (!strcmp(compname, "Analyzer_11")) return (void *) &(_Analyzer_11_var._parameters);
  if (!strcmp(compname, "Analyzer_12")) return (void *) &(_Analyzer_12_var._parameters);
  if (!strcmp(compname, "Analyzer_13")) return (void *) &(_Analyzer_13_var._parameters);
  if (!strcmp(compname, "Analyzer_14")) return (void *) &(_Analyzer_14_var._parameters);
  if (!strcmp(compname, "Analyzer_15")) return (void *) &(_Analyzer_15_var._parameters);
  if (!strcmp(compname, "Analyzer_16")) return (void *) &(_Analyzer_16_var._parameters);
  if (!strcmp(compname, "Analyzer_17")) return (void *) &(_Analyzer_17_var._parameters);
  if (!strcmp(compname, "Analyzer_18")) return (void *) &(_Analyzer_18_var._parameters);
  if (!strcmp(compname, "Analyzer_19")) return (void *) &(_Analyzer_19_var._parameters);
  if (!strcmp(compname, "Analyzer_20")) return (void *) &(_Analyzer_20_var._parameters);
  if (!strcmp(compname, "Analyzer_21")) return (void *) &(_Analyzer_21_var._parameters);
  if (!strcmp(compname, "Analyzer_22")) return (void *) &(_Analyzer_22_var._parameters);
  if (!strcmp(compname, "Analyzer_23")) return (void *) &(_Analyzer_23_var._parameters);
  if (!strcmp(compname, "Analyzer_24")) return (void *) &(_Analyzer_24_var._parameters);
  if (!strcmp(compname, "Analyzer_25")) return (void *) &(_Analyzer_25_var._parameters);
  if (!strcmp(compname, "Analyzer_26")) return (void *) &(_Analyzer_26_var._parameters);
  if (!strcmp(compname, "Analyzer_27")) return (void *) &(_Analyzer_27_var._parameters);
  if (!strcmp(compname, "Analyzer_28")) return (void *) &(_Analyzer_28_var._parameters);
  if (!strcmp(compname, "Analyzer_29")) return (void *) &(_Analyzer_29_var._parameters);
  if (!strcmp(compname, "Analyzer_30")) return (void *) &(_Analyzer_30_var._parameters);
  if (!strcmp(compname, "Analyzer_31")) return (void *) &(_Analyzer_31_var._parameters);
  if (!strcmp(compname, "Analyzer_32")) return (void *) &(_Analyzer_32_var._parameters);
  if (!strcmp(compname, "Analyzer_33")) return (void *) &(_Analyzer_33_var._parameters);
  if (!strcmp(compname, "Analyzer_34")) return (void *) &(_Analyzer_34_var._parameters);
  if (!strcmp(compname, "Analyzer_35")) return (void *) &(_Analyzer_35_var._parameters);
  if (!strcmp(compname, "Analyzer_36")) return (void *) &(_Analyzer_36_var._parameters);
  if (!strcmp(compname, "Analyzer_37")) return (void *) &(_Analyzer_37_var._parameters);
  if (!strcmp(compname, "Analyzer_38")) return (void *) &(_Analyzer_38_var._parameters);
  if (!strcmp(compname, "Analyzer_39")) return (void *) &(_Analyzer_39_var._parameters);
  if (!strcmp(compname, "Analyzer_40")) return (void *) &(_Analyzer_40_var._parameters);
  if (!strcmp(compname, "Analyzer_41")) return (void *) &(_Analyzer_41_var._parameters);
  if (!strcmp(compname, "Analyzer_42")) return (void *) &(_Analyzer_42_var._parameters);
  if (!strcmp(compname, "Analyzer_43")) return (void *) &(_Analyzer_43_var._parameters);
  if (!strcmp(compname, "Analyzer_44")) return (void *) &(_Analyzer_44_var._parameters);
  if (!strcmp(compname, "Analyzer_45")) return (void *) &(_Analyzer_45_var._parameters);
  if (!strcmp(compname, "rotator")) return (void *) &(_rotator_var._parameters);
  if (!strcmp(compname, "detector_pos")) return (void *) &(_detector_pos_var._parameters);
  if (!strcmp(compname, "detector")) return (void *) &(_detector_var._parameters);
  if (!strcmp(compname, "psd_monitor_end")) return (void *) &(_psd_monitor_end_var._parameters);
  if (!strcmp(compname, "E_PSD_mon_end")) return (void *) &(_E_PSD_mon_end_var._parameters);
  return 0;
}

void* _get_particle_var(char *token, _class_particle *p)
/* enables setpars based use of GET_PARTICLE_DVAR macro and similar */
{
  return 0;
}

int _getcomp_index(char* compname)
/* Enables retrieving the component position & rotation when the index is not known.
 * Component indexing into MACROS, e.g., POS_A_COMP_INDEX, are 1-based! */
{
  if (!strcmp(compname, "origin")) return 1;
  if (!strcmp(compname, "source")) return 2;
  if (!strcmp(compname, "psd_test")) return 3;
  if (!strcmp(compname, "e_monitor")) return 4;
  if (!strcmp(compname, "Analyzer_0")) return 5;
  if (!strcmp(compname, "Analyzer_1")) return 6;
  if (!strcmp(compname, "Analyzer_2")) return 7;
  if (!strcmp(compname, "Analyzer_3")) return 8;
  if (!strcmp(compname, "Analyzer_4")) return 9;
  if (!strcmp(compname, "Analyzer_5")) return 10;
  if (!strcmp(compname, "Analyzer_6")) return 11;
  if (!strcmp(compname, "Analyzer_7")) return 12;
  if (!strcmp(compname, "Analyzer_8")) return 13;
  if (!strcmp(compname, "Analyzer_9")) return 14;
  if (!strcmp(compname, "Analyzer_10")) return 15;
  if (!strcmp(compname, "Analyzer_11")) return 16;
  if (!strcmp(compname, "Analyzer_12")) return 17;
  if (!strcmp(compname, "Analyzer_13")) return 18;
  if (!strcmp(compname, "Analyzer_14")) return 19;
  if (!strcmp(compname, "Analyzer_15")) return 20;
  if (!strcmp(compname, "Analyzer_16")) return 21;
  if (!strcmp(compname, "Analyzer_17")) return 22;
  if (!strcmp(compname, "Analyzer_18")) return 23;
  if (!strcmp(compname, "Analyzer_19")) return 24;
  if (!strcmp(compname, "Analyzer_20")) return 25;
  if (!strcmp(compname, "Analyzer_21")) return 26;
  if (!strcmp(compname, "Analyzer_22")) return 27;
  if (!strcmp(compname, "Analyzer_23")) return 28;
  if (!strcmp(compname, "Analyzer_24")) return 29;
  if (!strcmp(compname, "Analyzer_25")) return 30;
  if (!strcmp(compname, "Analyzer_26")) return 31;
  if (!strcmp(compname, "Analyzer_27")) return 32;
  if (!strcmp(compname, "Analyzer_28")) return 33;
  if (!strcmp(compname, "Analyzer_29")) return 34;
  if (!strcmp(compname, "Analyzer_30")) return 35;
  if (!strcmp(compname, "Analyzer_31")) return 36;
  if (!strcmp(compname, "Analyzer_32")) return 37;
  if (!strcmp(compname, "Analyzer_33")) return 38;
  if (!strcmp(compname, "Analyzer_34")) return 39;
  if (!strcmp(compname, "Analyzer_35")) return 40;
  if (!strcmp(compname, "Analyzer_36")) return 41;
  if (!strcmp(compname, "Analyzer_37")) return 42;
  if (!strcmp(compname, "Analyzer_38")) return 43;
  if (!strcmp(compname, "Analyzer_39")) return 44;
  if (!strcmp(compname, "Analyzer_40")) return 45;
  if (!strcmp(compname, "Analyzer_41")) return 46;
  if (!strcmp(compname, "Analyzer_42")) return 47;
  if (!strcmp(compname, "Analyzer_43")) return 48;
  if (!strcmp(compname, "Analyzer_44")) return 49;
  if (!strcmp(compname, "Analyzer_45")) return 50;
  if (!strcmp(compname, "rotator")) return 51;
  if (!strcmp(compname, "detector_pos")) return 52;
  if (!strcmp(compname, "detector")) return 53;
  if (!strcmp(compname, "psd_monitor_end")) return 54;
  if (!strcmp(compname, "E_PSD_mon_end")) return 55;
  return -1;
}

/* embedding file "metadata-r.c" */

/** --- Contents of  metadata-r.c ---------------------------------------------------------------------------------- */
// Created by Gregory Tucker, Data Management Software Centre, European Spallation Source ERIC on 07/07/23.
#ifndef MCCODE_NAME
#include "metadata-r.h"
#endif

char * metadata_table_key_component(char* key){
  if (strlen(key) == 0) return NULL;
  char sep[2] = ":\0"; // matches any number of repeated colons
  // look for the separator in the provided key; strtok is allowed to modify the string, so copy it
  char * tok = malloc((strlen(key) + 1) * sizeof(char));
  strcpy(tok, key);
  char * pch = strtok(tok, sep); // this *is* the component name (if provided) -- but we need to move the pointer
  char * comp = malloc((1 + strlen(pch)) * sizeof(char));
  strcpy(comp, pch);
  if (tok) free(tok);
  return comp;
}
char * metadata_table_key_literal(char * key){
  if (strlen(key) == 0) return NULL;
  char sep[3] = ":\0";
  char * tok = malloc((strlen(key) + 1 ) * sizeof(char));
  strcpy(tok, key);
  char * pch = strtok(tok, sep); // this *is* the component name (if provided)
  if (pch) pch = strtok(NULL, sep); // either NULL or the literal name
  char * name = NULL;
  if (pch) {
    name = malloc((1 + strlen(pch)) * sizeof(char));
    strcpy(name, pch);
  }
  if (tok) free(tok);
  return name;
}
int metadata_table_defined(int no, metadata_table_t * tab, char * key){
  if (strlen(key) == 0){
    return no;
  }
  char * comp = metadata_table_key_component(key);
  char * name = metadata_table_key_literal(key);
  // look through the table for the matching component and literal names
  int number = 0;
  for (int i=0; i<no; ++i){
    if (!strcmp(comp, tab[i].source)){
      if (name == NULL || !strcmp(name, tab[i].name)) ++number;
    }
  }
  if (comp) free(comp);
  if (name) free(name);
  return number;
}
char * metadata_table_type(int no, metadata_table_t * tab, char * key){
  if (strlen(key) == 0) {
    fprintf(stderr, "Unable to check type of non-existent key\n");
    exit(1);
  }
  char * comp = metadata_table_key_component(key);
  char * name = metadata_table_key_literal(key);
  if (name == NULL){
    fprintf(stderr, "Unable to check type of literal for component %s without its name\n", comp);
    free(comp);
    exit(1);
  }
  char * type = NULL;
  for (int i=0; i<no; ++i){
    if (!strcmp(comp, tab[i].source) && !strcmp(name, tab[i].name)) type = tab[i].type;
  }
  if (comp) free(comp);
  if (name) free(name);
  return type;
}

char * metadata_table_literal(int no, metadata_table_t * tab, char * key){
  if (strlen(key) == 0) {
    fprintf(stderr, "Unable to retrieve literal for non-existent key\n");
    exit(1);
  }
  char * comp = metadata_table_key_component(key);
  char * name = metadata_table_key_literal(key);
  if (name == NULL){
    fprintf(stderr, "Unable to retrieve literal for component %s without its name\n", comp);
    free(comp);
    exit(1);
  }
  char * type = NULL;
  for (int i=0; i<no; ++i){
    if (!strcmp(comp, tab[i].source) && !strcmp(name, tab[i].name)) type = tab[i].value;
  }
  if (comp) free(comp);
  if (name) free(name);
  return type;
}
void metadata_table_print_all_keys(int no, metadata_table_t * tab){
  for (int i=0; i<no; ++i){
    printf("%s::%s ", tab[i].source, tab[i].name);
  }
  printf("\n");
}
int metadata_table_print_all_components(int no, metadata_table_t * tab){
  int count = 0;
  char ** known = malloc(no * sizeof(char*));
  for (int i=0; i<no; ++i){
    int unknown = 1;
    for (int j=0; j<count; ++j) if (!strcmp(tab[i].source, known[j])) unknown = 0;
    if (unknown) known[count++] = tab[i].source;
  }
  size_t nchar = 0;
  for (int i=0; i<count; ++i) nchar += strlen(known[i]) + 1;
  char * line = malloc((nchar + 1) * sizeof(char));
  line[0] = '\0';
  for (int i=0; i<count; ++i) sprintf(line, "%s%s ", line, known[i]);
  line[strlen(line)] = '\0'; // eat the trailing space
  printf("%s\n", line);
  free(line);
  free(known);
  return count;
}
int metadata_table_print_component_keys(int no, metadata_table_t * tab, char * key){
  char * comp = metadata_table_key_component(key);
  char * name = metadata_table_key_literal(key);
  int count = 0;
  for (int i=0; i<no; ++i) if (!strcmp(tab[i].source, comp) && (name == NULL || !strcmp(tab[i].name, name))) {
    if (name == NULL) printf("%s ", tab[i].name);
    ++count;
  }
  if (name != NULL) printf("%d", count); // replace count by strlen(tab[i].value)?
  printf("\n");
  return count;
}
/* -------------------------------------------------------------------------------------Contents of  metadata-r.c --- */
/* End of file "metadata-r.c". */

/* embedding file "mccode_main.c" */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
  /*  double run_num = 0; */
  time_t  t;
  clock_t ct;

#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, instrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

  ct = clock();

  // device and host functional RNG seed
  struct timeval tm;
  gettimeofday(&tm, NULL);
  mcseed = (long) tm.tv_sec*1000000 + tm.tv_usec;
  mcstartdate = (long)tm.tv_sec;  /* set start date before parsing options and creating sim file */
  // init global _particle.randstate for random number use
  // during init(), finally() and display(). NOTE: during trace, a local
  // "_particle" variable is present and thus used instead.
  srandom(_hash(mcseed-1));

#ifdef USE_MPI
  /* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      instrument_name, instrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per noe */
  }
#endif /* USE_MPI */

#ifdef OPENACC
#ifdef USE_MPI
  int num_devices = acc_get_num_devices(acc_device_nvidia);
  if(num_devices>0){
    int my_device = mpi_node_rank % num_devices;
    acc_set_device_num( my_device, acc_device_nvidia );
    printf("Have found %d GPU devices on rank %d. Will use device %d.\n", num_devices, mpi_node_rank, my_device);
  }else{
    printf("There was an issue probing acc_get_num_devices, fallback to host\n");
    acc_set_device_type( acc_device_host );
  }
#endif
#endif

  /* *** parse options ******************************************************* */
  SIG_MESSAGE("[" __FILE__ "] main START");
  mcformat = getenv(FLAVOR_UPPER "_FORMAT") ?
             getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  instrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */


#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif


/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */


  // init executed by master/host
  siminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("[" __FILE__ "] main INITIALISE");
  init();


#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#ifdef USE_MPI
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

// MT specific init, note that per-ray init is empty
#if RNG_ALG == 2
  mt_srandom(mcseed);
#endif


// main raytrace work loop
#ifndef FUNNEL
  // legacy version
  raytrace_all(mcncount, mcseed);
#else
  MPI_MASTER(
  // "funneled" version in which propagation is more parallelizable
  printf("\nNOTE: CPU COMPONENT grammar activated:\n 1) \"FUNNEL\" raytrace algorithm enabled.\n 2) Any SPLIT's are dynamically allocated based on available buffer size. \n");
	     );
  raytrace_all_funnel(mcncount, mcseed);
#endif


#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif


  // save/finally executed by master node/thread/host
  finally();


#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */


  return 0;
} /* mccode_main */
/* End of file "mccode_main.c". */

/* end of generated C code WARP.c */
