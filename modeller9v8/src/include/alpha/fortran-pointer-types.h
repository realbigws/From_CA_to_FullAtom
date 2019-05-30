/** \file alpha/fortran-pointer-types.h
 *  \brief Dummy types to handle Fortran pointers.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 *
 *             Auto-generated! Any changes will be lost.
 */

#ifndef MOD_POINT_TYPES_H
#define MOD_POINT_TYPES_H

#define F_PRIM1_ARRAY_SIZE 9
#define F_PRIM2_ARRAY_SIZE 12
#define F_PRIM3_ARRAY_SIZE 15
#define F_PRIM4_ARRAY_SIZE 18
#define F_CHAR1_ARRAY_SIZE 9
#define F_CHAR2_ARRAY_SIZE 12
#define F_CHAR3_ARRAY_SIZE 15
#define F_DERV1_ARRAY_SIZE 9
#define F_DERV2_ARRAY_SIZE 12
#define F_DERV3_ARRAY_SIZE 15
#define F_DERV_ARRAY_SIZE 1

#ifdef __cplusplus
extern "C" {
#endif

/** Fortran string length */
typedef long int_str;

struct mod_int1_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_int2_array { void * dummy[F_PRIM2_ARRAY_SIZE]; };
struct mod_int3_array { void * dummy[F_PRIM3_ARRAY_SIZE]; };
struct mod_int641_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_int642_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_int643_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_float1_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_float2_array { void * dummy[F_PRIM2_ARRAY_SIZE]; };
struct mod_float3_array { void * dummy[F_PRIM3_ARRAY_SIZE]; };
struct mod_float4_array { void * dummy[F_PRIM4_ARRAY_SIZE]; };
struct mod_double1_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_double2_array { void * dummy[F_PRIM2_ARRAY_SIZE]; };
struct mod_double3_array { void * dummy[F_PRIM3_ARRAY_SIZE]; };
struct mod_bool1_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_bool2_array { void * dummy[F_PRIM2_ARRAY_SIZE]; };
struct mod_bool3_array { void * dummy[F_PRIM3_ARRAY_SIZE]; };
struct mod_uchar1_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_uchar2_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_uchar3_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_short1_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_short2_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_short3_array { void * dummy[F_PRIM1_ARRAY_SIZE]; };
struct mod_char1_array { void * dummy[F_CHAR1_ARRAY_SIZE]; };
struct mod_char2_array { void * dummy[F_CHAR2_ARRAY_SIZE]; };
struct mod_char3_array { void * dummy[F_CHAR3_ARRAY_SIZE]; };
struct mod_derv1_array { void * dummy[F_DERV1_ARRAY_SIZE]; };
struct mod_derv2_array { void * dummy[F_DERV2_ARRAY_SIZE]; };
struct mod_derv3_array { void * dummy[F_DERV3_ARRAY_SIZE]; };
struct mod_derv_pt { void * dummy[F_DERV_ARRAY_SIZE]; };

#ifdef __cplusplus
}
#endif
#endif  /* MOD_POINT_TYPES_H */
