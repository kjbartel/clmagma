/*
    -- MAGMA (version 1.3.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2014

       @author Mark Gates
       @precisions normal z -> s d c

       auto-converted from zlanhe.cu

*/
#include "kernels_header.h"
#include "zlanhe.h"


/* ====================================================================== */
/* inf-norm */

/* Computes row sums dwork[i] = sum( abs( A(i,:) )), i=0:n-1, for || A ||_inf,
 * where n is any size and A is stored lower.
 * Has ceil( n / inf_bs ) blocks of (inf_bs x 4) threads each (inf_bs=32).
 * z precision uses > 16 KB shared memory, so requires Fermi (arch >= 200). */
__kernel void
zlanhe_inf_kernel_generic_lower(
    int n, const __global magmaDoubleComplex* A, unsigned long A_offset, int lda, 
    __global double *dwork, unsigned long dwork_offset,
    int n_full_block, int n_mod_bs )
{
    A += A_offset;
    dwork += dwork_offset;

    int tx = get_local_id(0);
    int ty = get_local_id(1);
    
    int diag = get_group_id(0)*inf_bs;
    int ind  = get_group_id(0)*inf_bs + tx;
    
    double res = 0.;
    
    __local magmaDoubleComplex la[inf_bs][inf_bs+1];
    
    if ( get_group_id(0) < n_full_block ) {
        // ------------------------------
        // All full block rows
        A += ind;
        A += ty * lda;
        
        // ----------
        // loop over all blocks left of the diagonal block
        for(int i=0; i < diag; i += inf_bs ) {
            // 32x4 threads cooperatively load 32x32 block
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                la[tx][ty+j] = A[j*lda];
            }
            A += lda*inf_bs;
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // compute 4 partial sums of each row, i.e.,
            // for ty=0:  res = sum( la[tx, 0: 7] )
            // for ty=1:  res = sum( la[tx, 8:15] )
            // for ty=2:  res = sum( la[tx,16:23] )
            // for ty=3:  res = sum( la[tx,24:31] )
            #pragma unroll 8             
            for(int j=ty*8; j < ty*8 + 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // load diagonal block
        #pragma unroll 8
        for(int j=0; j < inf_bs; j += 4) {
            la[tx][ty+j] = A[j*lda];
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // copy lower triangle to upper triangle, and
        // make diagonal real (zero imaginary part)
        #pragma unroll 8
        for(int i=ty*8; i < ty*8 + 8; i++) {
            if ( i < tx ) {
                la[i][tx] = la[tx][i];
            }
            #if defined(PRECISION_z) || defined(PRECISION_c)
            else if ( i == tx ) {
                la[i][i] = MAGMA_Z_MAKE( MAGMA_Z_REAL( la[i][i] ), 0 );
            }
            #endif
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // partial row sums
        #pragma unroll 8
        for(int j=ty*8; j < ty*8 + 8; j++) {
            res += MAGMA_Z_ABS( la[tx][j] );
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // ----------
        // loop over all 32x32 blocks below diagonal block
        A += inf_bs;
        for(int i=diag + inf_bs; i < n - n_mod_bs; i += inf_bs ) {
            // load block (transposed)
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                la[ty+j][tx] = A[j*lda];
            }
            A += inf_bs;
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // partial row sums
            #pragma unroll 8
            for(int j=ty*8; j < ty*8 + 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // last partial block, which is (n_mod_bs by inf_bs)
        if ( n_mod_bs > 0 ) {
            // load block (transposed), with zeros for rows outside matrix
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                if ( tx < n_mod_bs ) {
                    la[ty+j][tx] = A[j*lda];
                }
                else {
                    la[ty+j][tx] = MAGMA_Z_ZERO;
                }
            }
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // partial row sums
            #pragma unroll 8
            for(int j=ty*8; j < ty*8 + 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // 32x4 threads store partial sums into shared memory
        la[tx][ty] = MAGMA_Z_MAKE( res, 0. );
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // first column of 32x1 threads computes final sum of each row
        if ( ty == 0 ) {
            res = res
                + MAGMA_Z_REAL( la[tx][1] )
                + MAGMA_Z_REAL( la[tx][2] )
                + MAGMA_Z_REAL( la[tx][3] );
            dwork[ind] = res;
        }
    }
    else {
        // ------------------------------
        // Last, partial block row
        // Threads past end of matrix (i.e., ind >= n) are redundantly assigned
        // the last row (n-1). At the end, those results are ignored -- only
        // results for ind < n are saved into dwork.
        if ( tx < n_mod_bs ) {
            A += ind;
        }
        else {
            A += (get_group_id(0)*inf_bs + n_mod_bs - 1);  // redundantly do last row
        }
        A += ty * lda;
        
        // ----------
        // loop over all blocks left of the diagonal block
        // each is (n_mod_bs by inf_bs)
        for(int i=0; i < diag; i += inf_bs ) {
            // load block
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                la[tx][ty+j] = A[j*lda];
            }
            A += lda*inf_bs;
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // partial row sums
            #pragma unroll 8
            for(int j=0; j < 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j+ty*8] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // partial diagonal block
        if ( ty == 0 && tx < n_mod_bs ) {
            // sum rows left of diagonal
            for(int j=0; j < tx; j++) {
                res += MAGMA_Z_ABS( *A );
                A += lda;
            }
            // sum diagonal (ignoring imaginary part)
            res += MAGMA_D_ABS( MAGMA_Z_REAL( *A ));
            A += 1;
            // sum column below diagonal
            for(int j=tx+1; j < n_mod_bs; j++) {
                res += MAGMA_Z_ABS( *A );
                A += 1;
            }
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // ----------
        // 32x4 threads store partial sums into shared memory
        la[tx][ty]= MAGMA_Z_MAKE( res, 0. );
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // first column of 32x1 threads computes final sum of each row
        // rows outside matrix are ignored
        if ( ty == 0 && tx < n_mod_bs ) {
            res = res
                + MAGMA_Z_REAL( la[tx][1] )
                + MAGMA_Z_REAL( la[tx][2] )
                + MAGMA_Z_REAL( la[tx][3] );
            dwork[ind] = res;
        }
    }
}



/* Computes row sums dwork[i] = sum( abs( A(i,:) )), i=0:n-1, for || A ||_inf,
 * where n is any size and A is stored upper.
 * Has ceil( n / inf_bs ) blocks of (inf_bs x 4) threads each (inf_bs=32).
 * z precision uses > 16 KB shared memory, so requires Fermi (arch >= 200).
 * The upper implementation is similar to lower, but processes blocks
 * in the transposed order:
 * lower goes from left over to diagonal, then down to bottom;
 * upper goes from top  down to diagonal, then over to right.
 * Differences are noted with # in comments. */
__kernel void
zlanhe_inf_kernel_generic_upper(
    int n, const __global magmaDoubleComplex* A, unsigned long A_offset, int lda, 
    __global double *dwork, unsigned long dwork_offset,
    int n_full_block, int n_mod_bs )
{
    A += A_offset;
    dwork += dwork_offset;

    int tx = get_local_id(0);
    int ty = get_local_id(1);
    
    int diag = get_group_id(0)*inf_bs;
    int ind  = get_group_id(0)*inf_bs + tx;
    
    double res = 0.;
    
    __local magmaDoubleComplex la[inf_bs][inf_bs+1];
    
    if ( get_group_id(0) < n_full_block ) {
        // ------------------------------
        // All full block #columns
        A += get_group_id(0)*inf_bs*lda + tx;               //#
        A += ty * lda;
        
        // ----------
        // loop over all blocks #above the diagonal block
        for(int i=0; i < diag; i += inf_bs ) {
            // 32x4 threads cooperatively load 32x32 block (#transposed)
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                la[ty+j][tx] = A[j*lda];               //#
            }
            A += inf_bs;                               //#
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // compute 4 partial sums of each row, i.e.,
            // for ty=0:  res = sum( la[tx, 0: 7] )
            // for ty=1:  res = sum( la[tx, 8:15] )
            // for ty=2:  res = sum( la[tx,16:23] )
            // for ty=3:  res = sum( la[tx,24:31] )
            #pragma unroll 8             
            for(int j=ty*8; j < ty*8 + 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // load diagonal block
        #pragma unroll 8
        for(int j=0; j < inf_bs; j += 4) {
            la[tx][ty+j] = A[j*lda];
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // copy #upper triangle to #lower triangle, and
        // make diagonal real (zero imaginary part)
        #pragma unroll 8
        for(int i=ty*8; i < ty*8 + 8; i++) {
            if ( i > tx ) {                            //#
                la[i][tx] = la[tx][i];
            }
            #if defined(PRECISION_z) || defined(PRECISION_c)
            else if ( i == tx ) {
                la[i][i] = MAGMA_Z_MAKE( MAGMA_Z_REAL( la[i][i] ), 0 );
            }
            #endif
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // partial row sums
        #pragma unroll 8
        for(int j=ty*8; j < ty*8 + 8; j++) {
            res += MAGMA_Z_ABS( la[tx][j] );
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // ----------
        // loop over all 32x32 blocks #right of diagonal block
        A += inf_bs*lda;                               //#
        for(int i=diag + inf_bs; i < n - n_mod_bs; i += inf_bs ) {
            // load block (#non-transposed)
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                la[tx][ty+j] = A[j*lda];               //#
            }
            A += inf_bs*lda;                           //#
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // partial row sums
            #pragma unroll 8
            for(int j=ty*8; j < ty*8 + 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // last partial block, which is #(inf_bs by n_mod_bs)
        if ( n_mod_bs > 0 ) {
            // load block (#non-transposed), with zeros for #cols outside matrix
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                if ( ty+j < n_mod_bs ) {               //#
                    la[tx][ty+j] = A[j*lda];           //#
                }
                else {
                    la[tx][ty+j] = MAGMA_Z_ZERO;       //#
                }
            }
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // partial row sums
            #pragma unroll 8
            for(int j=ty*8; j < ty*8 + 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // 32x4 threads store partial sums into shared memory
        la[tx][ty] = MAGMA_Z_MAKE( res, 0. );
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // first column of 32x1 threads computes final sum of each row
        if ( ty == 0 ) {
            res = res
                + MAGMA_Z_REAL( la[tx][1] )
                + MAGMA_Z_REAL( la[tx][2] )
                + MAGMA_Z_REAL( la[tx][3] );
            dwork[ind] = res;
        }
    }
    else {
        // ------------------------------
        // Last, partial block #column
        // Instead of assigning threads ind >= n to the last row (n-1), as in Lower,
        // Upper simply adjusts loop bounds to avoid loading columns outside the matrix.
        // Again, at the end, those results are ignored -- only
        // results for ind < n are saved into dwork.
        A += get_group_id(0)*inf_bs*lda + tx;               //#
        A += ty * lda;
        
        // ----------
        // loop over all blocks #above the diagonal block
        // each is #(inf_bs by n_mod_bs)
        for(int i=0; i < diag; i += inf_bs ) {
            // load block (#transposed), #ignoring columns outside matrix
            #pragma unroll 8
            for(int j=0; j < inf_bs; j += 4) {
                if ( ty+j < n_mod_bs ) {
                    la[ty+j][tx] = A[j*lda];
                }
            }
            A += inf_bs;                               //#
            barrier( CLK_LOCAL_MEM_FENCE );
            
            // partial row sums
            #pragma unroll 8
            for(int j=0; j < 8; j++) {
                res += MAGMA_Z_ABS( la[tx][j+ty*8] );
            }
            barrier( CLK_LOCAL_MEM_FENCE );
        }
        
        // ----------
        // partial diagonal block
        if ( ty == 0 && tx < n_mod_bs ) {
            // #transpose pointer within diagonal block
            // #i.e., from A = A(tx,ty), transpose to A = A(ty,tx).
            A = A - tx - ty*lda + tx*lda + ty;
            
            // sum #column above diagonal
            for(int j=0; j < tx; j++) {
                res += MAGMA_Z_ABS( *A );
                A += 1;                                //#
            }
            // sum diagonal (ignoring imaginary part)
            res += MAGMA_D_ABS( MAGMA_Z_REAL( *A ));
            A += lda;                                  //#
            // sum #row right of diagonal
            for(int j=tx+1; j < n_mod_bs; j++) {
                res += MAGMA_Z_ABS( *A );
                A += lda;                              //#
            }
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // ----------
        // 32x4 threads store partial sums into shared memory
        la[tx][ty]= MAGMA_Z_MAKE( res, 0. );
        barrier( CLK_LOCAL_MEM_FENCE );
        
        // first column of 32x1 threads computes final sum of each row
        // rows outside matrix are ignored
        if ( ty == 0 && tx < n_mod_bs ) {
            res = res
                + MAGMA_Z_REAL( la[tx][1] )
                + MAGMA_Z_REAL( la[tx][2] )
                + MAGMA_Z_REAL( la[tx][3] );
            dwork[ind] = res;
        }
    }
}


/* ====================================================================== */
/* max-norm */

/* Computes dwork[i] = max( abs( A(i,0:i) )), i=0:n-1, for ||A||_max, where A is stored lower */
__kernel void
zlanhe_max_kernel_lower(
    int n, const __global magmaDoubleComplex* A, unsigned long A_offset, int lda, 
    __global double *dwork, unsigned long dwork_offset )
{
    A += A_offset;
    dwork += dwork_offset;

    int ind = get_group_id(0)*max_bs + get_local_id(0);
    double res = 0;

    if (ind < n) {
        A += ind;
        for(int j=0; j < ind; ++j) {
            res = fmax( res, MAGMA_Z_ABS( *A ));
            A += lda;
        }
        // diagonal element (ignoring imaginary part)
        res = fmax( res, MAGMA_D_ABS( MAGMA_Z_REAL( *A )));
        dwork[ind] = res;
    }
}


/* Computes dwork[i] = max( abs( A(i,0:i) )), i=0:n-1, for ||A||_max, where A is stored upper. */
__kernel void
zlanhe_max_kernel_upper(
    int n, const __global magmaDoubleComplex* A, unsigned long A_offset, int lda, 
    __global double *dwork, unsigned long dwork_offset )
{
    A += A_offset;
    dwork += dwork_offset;

    int ind = get_group_id(0)*max_bs + get_local_id(0);
    double res = 0;

    if (ind < n) {
        A += ind;
        A += (n-1)*lda;
        for(int j=n-1; j > ind; j--) {
            res = fmax( res, MAGMA_Z_ABS( *A ));
            A -= lda;
        }
        // diagonal element (ignoring imaginary part)
        res = fmax( res, MAGMA_D_ABS( MAGMA_Z_REAL( *A )));
        dwork[ind] = res;
    }
}
