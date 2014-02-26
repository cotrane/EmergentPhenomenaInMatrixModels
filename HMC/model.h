
/*
 * macro to define which model to simulate:
 *    MM8 is d=8 model
 *    MM3 is d=3 model
 */
#define MM8

#ifdef MM8
#define NUMMAT	8
#define LIEGROUP 3
#define	action_function	action8MM
#define deltaaction_function deltaX8MM
#define genmom gen_gaussmomcplx
#define compmom mom
#elif defined(MM3)
#define NUMMAT	3
#define LIEGROUP 2
#define	action_function	action3MMwrap
#define deltaaction_function deltaX3MMwrap
#define genmom gen_gaussmomcplx
#define compmom mom
#else
#define NUMMAT	8
#define LIEGROUP 3
#define	action_function	action8MM
#define deltaaction_function deltaX8MM
#define genmom gen_gaussmomcplx
#define compmom mom
#endif
