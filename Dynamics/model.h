
/*
 * change the define value below to change the model that is simulated
 * MM1 ... 1-dimensional Harmonic Oscillator
 * MM3 ... 3-d Yang-Mills-Myers model
 * MM8 ... 8-d Yang-Mills-Myers model
 */

#define MM8

#ifdef MM8
#define NUMMAT	8
#define LIEGROUP 3
#define	action_function	action8MM
#define deltaaction_function deltaX8MM
#define deltamom_function addmat
#define genmom gen_gaussmomcplx
#define compmom mom
#elif defined(MM3)
#define NUMMAT	3
#define LIEGROUP 2
#define	action_function	action3MMwrap
#define deltaaction_function deltaX3MMwrap
#define deltamom_function addmat
#define genmom gen_gaussmomcplx
#define compmom mom
#elif defined(MM1)
#define NUMMAT	1
#define LIEGROUP 1
#define	action_function	action1MMwrap
#define deltaaction_function deltaX1MMwrap
#define deltamom_function addmat1MM
#define genmom gen_gaussmomcplxtrace
#define compmom mom
#else
#define NUMMAT	8
#define LIEGROUP 3
#define	action_function	action8MM
#define deltaaction_function deltaX8MM
#define deltamom_function addmat
#define genmom gen_gaussmomcplx
#define compmom mom
#endif
