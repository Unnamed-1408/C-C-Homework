#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <immintrin.h>
#include <cstdlib>
#include "/usr/lib/gcc/x86_64-redhat-linux/10/include/omp.h"

using namespace std;

struct Matrix{
    size_t row = 0; //行
    size_t column = 0; //列
    float **Amn;
};
void Test();
void output_Matrix(Matrix a);
Matrix* Create_Matrix_One(size_t row,size_t column);
Matrix* Create_Matrix(Matrix* a, Matrix* b);
inline void Matrix_Multiplication(Matrix* a, Matrix* b, Matrix* output);
static inline void Matrix_Multiplication_order(Matrix* a, Matrix* b, Matrix* output);
void Matrix_Multiplication_avx(Matrix *A, Matrix *B, Matrix *C);
inline double sse_inner(const float* a, const float* b, unsigned int size);
void dgemm_avx_unroll_blk_omp(size_t n, float *A, float *B, float *C);
void Cache_Blocking(Matrix* a, Matrix* b, Matrix* output);
void matrix_mult_wiki_block(const float*A , const float* B, float* C,const int N, const int M, const int K) ;
void mul_matrix(Matrix *result, Matrix *mat1, Matrix *mat2);
static inline void Matrix_Multiplication_order_try(Matrix* a, Matrix* b, Matrix* output);
//inline void Matrix_Multiplication_sse(Matrix* a, Matrix* b, Matrix* output);
//void AVX_Test(Matrix* a, Matrix* b, Matrix* output);
//inline void Matrix_Multiplication_order_UNROLL(Matrix* a, Matrix* b, Matrix* output);



int main()
{
    string in = "";
    while(in != "q" && in != "quit"){
        Test();
        cout << "Enter q/quit to quit; Input any other to continue" <<endl;
        cin >> in;
    }
}

void Test(){
/*
    Matrix* a = Create_Matrix_One(10,10);
    Matrix* b = Create_Matrix_One(4000,4000);
    ifstream inF;
    inF.open("martix_a.dat", std::ifstream::binary);
    inF.read(reinterpret_cast<char*>(a->Amn), sizeof(float)*(10*10));
    inF.close();

    ifstream inFa;
    inFa.open("martix_b.dat", std::ifstream::binary);
    inFa.read(reinterpret_cast<char*>(b->Amn), sizeof(float)*(4000*4000));
    inFa.close();
    */
    int a_row = 0,a_column = 0;
    int b_row = 0,b_column = 0;
    cout << "Please enter the first matrix's row and column:" << endl;
    cout << "row:";
    cin >> a_row ;
    cout << "column:";
    cin >> a_column;
    cout << "The second one:" << endl;
    cout << "row:";
    cin >> b_row;
    cout << "column:";
    cin >> b_column;
    Matrix* a = Create_Matrix_One(a_row,a_column);
    Matrix* b = Create_Matrix_One(b_row,b_column);
    if(a_column != b_row){
        cout << "Check your input! These two matrix cannot multiply!" << endl;
        delete a,b;
        return;
    }

    cout << "Random or not? No random will put all the elements as 1 or input by yourself.Y/N" << endl;
    string adsa = "";
    string mm = "";
    cin >> adsa;
    bool randm = false;
    bool input = false;
    if(adsa == "Y"){
        randm = true;
        cout << "Random is setted ! Creating Matrix..." << endl;
    }
    else{
        cout << "Random is not setted ! Do you need to input the matrix by yourself?" << endl;
        cin >> mm;
        if(mm == "Y"){
            input = true;
            cout << "You choose to input yourself,The frist Matrix first" << endl;
        }
        else
            cout << "Creating Matrix..." << endl;
    }

    for(size_t i = 0; i<a->row;i++){
        for(size_t j = 0; j<a->column;j++){
            if(randm == false){
                if(input ==  false)
                    a->Amn[i][j] = 1;
                else{
                    cin >> a->Amn[i][j];
                }
            }
            else
                a->Amn[i][j] = (rand() / float(RAND_MAX));
        }
    }
    cout << "Making the second..." << endl;
    for(size_t i = 0; i<b->row;i++){
        for(size_t j = 0; j<b->column;j++){
            if(randm == false)
                if(input == false)
                    b->Amn[i][j] = 1;
                else{
                    cin >> b->Amn[i][j];
                }
            else
                b->Amn[i][j] = (rand() / float(RAND_MAX))*__FLT_MAX__;
        }
    }

    Matrix *output = Create_Matrix(a,b);

   auto start = std::chrono::steady_clock::now();
//code start here

//    Matrix_Multiplication(a,b,output);

   Matrix_Multiplication_order(a,b,output);
//   Cache_Blocking(a,b,output);
//   mul_matrix(output,a,b);
//   matrix_mult_wiki_block(*a->Amn , *b->Amn, *output->Amn,4000, 4000, 4000) ;
//   dgemm_avx_unroll_blk_omp(2048, *a->Amn, *a->Amn, *output->Amn);
//   AVX_Test(a,b,output);
//    Matrix_Multiplication_avx(a,b,output);

//code end here


    auto end = std::chrono::steady_clock::now();
    cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << endl;
    cout << "Do you want to output the result?" << endl;
    string endb = "";
    cin >> endb;
    if(endb == "Y")
        output_Matrix(*output);
}


inline void Matrix_Multiplication(Matrix* a, Matrix* b, Matrix* output){
#pragma omp parallel for
    for(size_t i = 0; i < output->row-8; i++){
        for(size_t j =0; j < output->column; j++){
            for(size_t k=0; k < a->column; k++){
                output->Amn[i][j] += a->Amn[i][k] * b->Amn[k][j];
            }
        }
    }
}


//inline void Matrix_Multiplication(Matrix* a, Matrix* b, Matrix* output){
//#pragma omp parallel for
//    for(size_t i = 0; i < output->row-8; i++){
//        for(size_t j =0; j < output->column; j++){
//            for(size_t k=0; k < a->column; k+=8){
//                __m256 veca = _mm256_set_ps(a->Amn[i][k],a->Amn[i][k+1],a->Amn[i][k+2],a->Amn[i][k+3],a->Amn[i][k+4],a->Amn[i][k+6],a->Amn[i][k+7],a->Amn[i][k+8]);
//                __m256 vecb = _mm256_set_ps(b->Amn[k][j],a->Amn[k+1][j],a->Amn[k+2][j],a->Amn[k+3][j],a->Amn[k+4][j],a->Amn[k+5][j],a->Amn[k+6][j],a->Amn[k+7][j]);
//                __m256 vecans = _mm256_mul_ps(veca,vecb);
//                float ans[8];
//                _mm256_storeu_ps(ans,vecans);
//                output->Amn[i][j] += ans[0]+ans[1]+ans[2]+ans[3]+ans[4]+ans[5]+ans[6]+ans[7];
//            }
//        }
//    }
//}


static inline void Matrix_Multiplication_order(Matrix* a, Matrix* b, Matrix* output){
#pragma omp parallel for
    for(size_t i = 0; i < output->row / 4 *4; i += 4){
        for(size_t k =0; k < a->column; k++){
            register float temp[4];
            temp[0] = a->Amn[i][k];
            temp[1] = a->Amn[i+1][k];
            temp[2] = a->Amn[i+2][k];
            temp[3] = a->Amn[i+3][k];
            register auto ac = &temp[0];
            for(size_t j=0; j < output->column; j++){
                auto cd = b->Amn[k][j];
                output->Amn[i][j] += *(ac) * cd;
                output->Amn[i+1][j] += *(ac+1) * cd;
                output->Amn[i+2][j] += *(ac+2) * cd;
                output->Amn[i+3][j] += *(ac+3) * cd;
            }
        }
    }
    for(size_t i = 0; i < output->row -output->row / 4 * 4; i++){
        for(size_t k =0; k < a->column; k++){
            float temp = a->Amn[output->row/4 *4 + i][k];
            for(size_t j=0; j < output->column; j++){
                auto cd = b->Amn[k][j];
                output->Amn[output->row /4 *4 + i][j] += temp * cd;
            }
        }
    }
}

#define BLOCK_SIZE 32
static inline void Matrix_Multiplication_order_try(Matrix* a, Matrix* b, Matrix* output){
#pragma omp parallel for
    for(size_t i = 0; i < output->row-4; i += BLOCK_SIZE){
        for(size_t k =0; k < a->column; k += BLOCK_SIZE){
            for(size_t j=0; j < output->column; j += BLOCK_SIZE){
                for(size_t ii = i; ii < i+BLOCK_SIZE; ii += 4){
                    for(size_t kk =k; kk < k+BLOCK_SIZE; kk ++){
                        register float temp[4];
                        temp[0] = a->Amn[ii][kk];
                        temp[1] = a->Amn[ii+1][kk];
                        temp[2] = a->Amn[ii+2][kk];
                        temp[3] = a->Amn[ii+3][kk];
                        register auto ac = &temp[0];
                        for(size_t jj=j; jj < j+BLOCK_SIZE; jj ++){
                            auto cd = b->Amn[kk][jj];
                            output->Amn[ii][jj] += *(ac) * cd;
                            output->Amn[ii+1][jj] += *(ac+1) * cd;
                            output->Amn[ii+2][jj] += *(ac+2) * cd;
                            output->Amn[ii+3][jj] += *(ac+3) * cd;

                        }
                    }
                }
            }
        }
    }
}






static inline void Matrix_Multiplication_order_blocking(Matrix* a, Matrix* b, Matrix* output , size_t si, size_t sj, size_t sk){
    for(size_t i = si; i < si+BLOCK_SIZE; i +=4){
        for(size_t k = sk; k < sk+BLOCK_SIZE; k++){
            register float temp[4];
            temp[0] = a->Amn[i][k];
            temp[1] = a->Amn[i+1][k];
            temp[2] = a->Amn[i+2][k];
            temp[3] = a->Amn[i+3][k];
            register auto ac = &temp[0];
            for(size_t j=sj; j < sj+BLOCK_SIZE; j++){
                auto cd = b->Amn[k][j];
                output->Amn[i][j] += *(ac) * cd;
                output->Amn[i+1][j] += *(ac+1) * cd;
                output->Amn[i+2][j] += *(ac+2) * cd;
                output->Amn[i+3][j] += *(ac+3) * cd;
            }
        }
    }
}



void Cache_Blocking(Matrix* a, Matrix* b, Matrix* output){
#pragma omp parallel for
    for(size_t i = 0; i < output->row; i += BLOCK_SIZE){
        for(size_t k =0; k < a->column; k += BLOCK_SIZE){
            for(size_t j=0; j < output->column; j += BLOCK_SIZE){
                Matrix_Multiplication_order_blocking(a,b,output,i,j,k);
            }
        }
    }
}


void mul_matrix(Matrix *result, Matrix *mat1, Matrix *mat2){
    int I = mat1->row;
    int J = mat2->column;
    int K = mat2->row;
    #pragma omp parallel for
    for(int i = 0; i < I; i++){
        for(int k = 0; k < K; k++){
            __m256 vA = _mm256_set1_ps((*mat1->Amn)[i * K + k]);
            for(int j = 0; j < J / 8 * 8; j += 8){
                __m256 sum = _mm256_loadu_ps(*result->Amn + i * J + j);
                __m256 vB = _mm256_loadu_ps(*mat2->Amn + k * J + j);
                __m256 intermediate = _mm256_mul_ps(vA, vB);
                sum = _mm256_add_ps(sum, intermediate);
                _mm256_storeu_ps(*result->Amn + i * J + j, sum);
             }
             for(int x = J / 8 * 8; x < J; x++){
                 *result->Amn[i * J + x] += *mat1 -> Amn[i * K + k] * *mat2 -> Amn[k * J + x];
             }
         }
     }
}


void matrix_mult_wiki_block(const float*A , const float* B, float* C, const int N, const int M, const int K) {
   const int block_size = 8;  //I have tried several different block sizes
   for(int i=0; i<N; i++) {
       for(int j=0; j<K; j++) {
           C[K*i + j] = 0;
       }
    }
    for(int l2=0; l2<M; l2+=block_size) {
        for(int j2=0; j2<K; j2+=block_size) {
        #pragma omp parallel for
            for(int i=0; i<N; i++) {
                for(int l=l2; l<min(M, l2+block_size); l++) {
                    for(int j=j2; j<min(K, j2+block_size); j++) {
                        C[K*i + j] += A[M*i+l]*B[K*l+j];
                    }
                }
            }
        }
    }
}

#define UNROLL 4
#define BLOCKSIZE 32

static inline void do_block(int n, int si, int sj, int sk, float *A, float *B, float *C)
{
    for (int i = si; i < si + BLOCKSIZE; i += UNROLL*4) {
        for (int j = sj; j < sj + BLOCKSIZE; j++) {
            __m256 c[UNROLL];
            for (int x = 0; x < UNROLL; x++) {
                c[x] = _mm256_loadu_ps(C+i+x*4+j*n);
            }
            for (int k = sk; k < sk + BLOCKSIZE; k++) {
                __m256 b = _mm256_broadcast_ss(B+k+j*n);
                for (int x = 0; x < UNROLL; x++) {
                    c[x] = _mm256_add_ps(c[x],_mm256_mul_ps(_mm256_loadu_ps(A+n*k+x*4+i), b));
                }
            }

            for (int x = 0; x < UNROLL; x++) {
                _mm256_storeu_ps(C+i+x*4+j*n, c[x]);
            }
        }
    }
}




void dgemm_avx_unroll_blk_omp(size_t n, float *A, float *B, float *C)
{
#pragma omp parallel for
    for (int sj = 0; sj < n; sj += BLOCKSIZE) {
        for (int si = 0; si < n; si += BLOCKSIZE) {
            for (int sk = 0; sk < n; sk += BLOCKSIZE) {
                do_block(n, si, sj, sk, A, B, C);
            }
        }
    }
}


//inline void Matrix_Multiplication_order(Matrix* a, Matrix* b, Matrix* output){
//#pragma omp parallel for
//    for(size_t i = 0; i < output->row-8; i += 8){
//        register float ans[8];
//        for(size_t k =0; k < a->column; k++){
//            __m256 veca = _mm256_set_ps(a->Amn[i][k],a->Amn[i+1][k],a->Amn[i+2][k],a->Amn[i+3][k],a->Amn[i+4][k],a->Amn[i+5][k],a->Amn[i+6][k],a->Amn[i+7][k]);
//            for(size_t j=0; j < output->column; j++){
//                __m256 vecb = _mm256_set1_ps(b->Amn[k][j]);
//                __m256 vecans = _mm256_mul_ps(veca,vecb);

//                _mm256_storeu_ps(ans,vecans);
//                output->Amn[i][j] += ans[0];
//                output->Amn[i+1][j] += ans[1];
//                output->Amn[i+2][j] += ans[2];
//                output->Amn[i+3][j] += ans[3];
//                output->Amn[i+4][j] += ans[4];
//                output->Amn[i+5][j] += ans[5];
//                output->Amn[i+6][j] += ans[6];
//                output->Amn[i+7][j] += ans[7];
//            }
//        }
//    }
//}
/*
inline void Matrix_Multiplication_sse(Matrix* a, Matrix* b, Matrix* output){
#pragma omp parallel for
    for(size_t i = 0; i < output->row; i = i+4){
            for(size_t j=0; j < output->column; j++){
                output->Amn[i][j] = sse_inner(a->Amn[i],b->Amn[j],4000);
                output->Amn[i+1][j] = sse_inner(a->Amn[i+1],b->Amn[j],4000);
                output->Amn[i+2][j] = sse_inner(a->Amn[i+2],b->Amn[j],4000);
                output->Amn[i+3][j] = sse_inner(a->Amn[i+3],b->Amn[j],4000);
            }
    }
}


inline double sse_inner(const float* a, const float* b, unsigned int size)
{
        float z = 0.0f;
        double fres = 0;
        float ftmp[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
        __m128 mres;

        if ((size / 4) != 0) {
                mres = _mm_load_ss(&z);
                for (unsigned int i = 0; i < size / 4; i++)
                        mres = _mm_add_ps(mres, _mm_mul_ps(_mm_loadu_ps(&a[4*i]),_mm_loadu_ps(&b[4*i])));


                __m128 mv1 = _mm_movelh_ps(mres, mres);
                __m128 mv2 = _mm_movehl_ps(mres, mres);     //c,d,c,d
                mres = _mm_add_ps(mv1, mv2);                //res[0],res[1]

                _mm_store_ps(ftmp, mres);

                fres = (double)ftmp[0] + (double)ftmp[1];
        }

        if ((size % 4) != 0) {
                for (unsigned int i = size - size % 4; i < size; i++)
                        fres += (double)a[i] * (double)b[i];
        }

        return fres;
}
*/
/*

#define UNROLL 4
#define BLOCKSIZE 32

inline void do_block(size_t si,size_t sj,size_t sk,Matrix* a, Matrix* b, Matrix* output){

    for (size_t i = si; i < si + BLOCKSIZE; i += UNROLL) {
        for (size_t j = sj; j < sj + BLOCKSIZE; j++) {
            for (size_t k = sk; k < sk + BLOCKSIZE; k++) {
                for (int x = 0; x < UNROLL; x++) {
                    output->Amn[i][j] +=a->Amn[i][k] * b->Amn[k][j];
                }
            }
        }
    }
}


inline void Matrix_Multiplication_order_UNROLL(Matrix* a, Matrix* b, Matrix* output){
    size_t i = 0;
#pragma omp parallel for
    for(i = 0; i < output->row-BLOCKSIZE; i += BLOCKSIZE){
        size_t k =0;
        for(k =0; k < a->column-BLOCKSIZE; k += BLOCKSIZE){
            size_t j = 0;
            for(j=0; j < output->column-BLOCKSIZE; j += BLOCKSIZE){
                do_block(i,j,k,a,b,output);
            }
            while(j != output->column){
                output->Amn[i][j] += a->Amn[i][k] * b->Amn[k][j];
                j++;
            }
        }
        while(k != a->column-BLOCKSIZE){
            size_t j = 0;
            for(j=0; j < output->column-BLOCKSIZE; j += BLOCKSIZE){
                do_block(i,j,k,a,b,output);
            }
            while(j != output->column){
                output->Amn[i][j] += a->Amn[i][k] * b->Amn[k][j];
                j++;
            }
        }
    }

    for(; i < output->row; i++){
        for(size_t k =0; k < a->column; k++){
            float temp = a->Amn[i][k];
            for(size_t j=0; j < output->column; j++){
                output->Amn[i][j] += temp * b->Amn[k][j];
            }
        }
    }


}

*/


void Matrix_Multiplication_avx(Matrix *A, Matrix *B, Matrix *C)
{
    size_t i = 0;
#pragma omp parallel for
    for(i = 0; i < C->row - 8; i += 8){
        for(size_t j =0; j < C->column; j++){
            __m256 c0 = _mm256_loadu_ps(&C->Amn[i][j]);
            for(size_t k=0; k < A->column; k++){
                c0 = _mm256_add_ps(c0,_mm256_mul_ps(_mm256_loadu_ps(&A->Amn[i][k]),_mm256_broadcast_ss(&B->Amn[k][j])));
            }
            _mm256_storeu_ps(&C->Amn[i][j],c0);
        }
    }

    if(i != C->row){
        for(; i < C->row; i++){
            for(size_t k =0; k < A->column; k++){
                float temp = A->Amn[i][k];
                for(size_t j=0; j < C->column; j++){
                    C->Amn[i][j] += temp * B->Amn[k][j];
                }
            }
        }
    }
}

/*
#define UNROLL 4
#define BLOCKSIZE 32

static inline void do_block(int si, int sj, int sk,	Matrix *A, Matrix *B, Matrix *C)
{
    for (int i = si; i < si + BLOCKSIZE; i += UNROLL*8) {
        for (int j = sj; j < sj + BLOCKSIZE; j++) {
            __m256 c[UNROLL];
            for (int x = 0; x < UNROLL; x++) {
//				c[x] = _mm256_loadu_ps(C+i+x*4+j*n);
            }
            for (int k = sk; k < sk + BLOCKSIZE; k++) {
//				__m256 b = _mm256_broadcast_ss(B+k+j*n);
                for (int x = 0; x < UNROLL; x++) {
//                    c[x] = _mm256_add_ps(c[x],_mm256_mul_ps(_mm256_loadu_ps(A+n*k+x*4+i), b));
                }
            }

            for (int x = 0; x < UNROLL; x++) {
//				_mm256_store_ps(C+i+x*4+j*n, c[x]);
            }
        }
    }
}

void dgemm_avx_unroll_blk_omp(Matrix *A, Matrix *B, Matrix *C)
{
#pragma omp parallel for
    for (size_t sj = 0; sj < C->column; sj += BLOCKSIZE) {
        for (size_t si = 0; si < C->row; si += BLOCKSIZE) {
            for (size_t sk = 0; sk < A->column; sk += BLOCKSIZE) {
                do_block(si, sj, sk, A, B, C);
            }
        }
    }
}
*/




/*
#define BLOCK 16
void AVX_Test(Matrix* a, Matrix* b, Matrix* output){
    size_t c0 = 0;
    size_t aWidth256(a->column / 8);
    size_t aWidthWarpFloor(aWidth256 / BLOCK * BLOCK);
    __m256* aData((__m256*)b->Amn);//乘法的右矩阵
    __m256* bData((__m256*)output->Amn);//乘法结果

//#pragma omp parallel for
    for(c0 = 0;c0 < a->row-a->row%2; c0 += 2){
        size_t c1 = 0;
        for(; c1 < aWidthWarpFloor; c1 +=BLOCK){
            __m256 ans0[BLOCK] = {0};
            __m256 ans1[BLOCK] = {0};
            for(size_t c2 = 0; c2 < a->column; c2++){
                __m256 tp0 = _mm256_set1_ps(a->Amn[c0][c2]);
                __m256 tp1 = _mm256_set1_ps(a->Amn[c0+1][c2]);

                for(size_t c3 = 0; c3 < BLOCK; c3++){
                    __m256 b = aData[aWidth256 * c2 + c1 + c3];
                    ans0[c3] = _mm256_fmadd_ps(tp0,b,ans0[c3]);
                    ans1[c3] = _mm256_fmadd_ps(tp1, b, ans1[c3]);
                }
            }

            for(size_t c3 = 0; c3 < BLOCK; c3++){
                bData[c0 * aWidth256 + c1 +c3] = ans0[c3];
                bData[(c0 + 1) * aWidth256 + c1 +c3] = ans1[c3];
            }
        }

        if (c1 < aWidth256){
            __m256 ans0[BLOCK] = {0};
            __m256 ans1[BLOCK] = {0};
            for(size_t c2 = 0; c2 < a->column; c2++){
                __m256 tp0 = _mm256_set1_ps(a->Amn[c0][c2]);
                __m256 tp1 = _mm256_set1_ps(a->Amn[c0+1][c2]);

                for(size_t c3 = 0; c3 < aWidth256-aWidthWarpFloor; c3++){
                    __m256 b = aData[aWidth256 * c2 + c1 + c3];
                    ans0[c3] = _mm256_fmadd_ps(tp0,b,ans0[c3]);
                    ans1[c3] = _mm256_fmadd_ps(tp1, b, ans1[c3]);
                }
            }

            for(size_t c3 = 0; c3 < aWidth256-aWidthWarpFloor; c3++){
                bData[c0 * aWidth256 + c1 +c3] = ans0[c3];
                bData[(c0 + 1) * aWidth256 + c1 +c3] = ans1[c3];
            }
        }
    }
}
*/

//Matrix* Create_Matrix(Matrix* a, Matrix* b){
//    size_t column = a->row;
//    size_t row = b->column;
//    Matrix* output = new Matrix;

//    output->row = row;
//    output->column = column;

//    output->Amn = new float*[row];
//    for(size_t i=0; i<row; ++i)
//    {
//       output->Amn[i]=new float[column]();
//    }

//    return output;
//}




Matrix* Create_Matrix(Matrix* a, Matrix* b){
    size_t row = a->row;
    size_t column = b->column;
    Matrix* output = new Matrix;

    output->row = row;
    output->column = column;

    output->Amn = new float*[row];
    for(size_t i=0; i<row; ++i)
    {
       output->Amn[i]=new float[column]();
    }

    return output;
}

Matrix* Create_Matrix_One(const size_t row,const size_t column){
    Matrix* output = new Matrix;

    output->row = row;
    output->column = column;
    output->Amn = new float*[row];
    for(size_t i=0; i<row; ++i)
    {
       output->Amn[i]=new float[column];
    }
    return output;
}


void output_Matrix(Matrix a){
    size_t temp_a = 0;
    size_t temp_b = 0;
    while(temp_a != a.row){
        temp_b = 0;
        while(temp_b != a.column){
            printf("%f ",a.Amn[temp_a][temp_b]);
//            printf("%f ",a.Amn[temp_a*a.row+temp_b]);
            temp_b++;
        }
        temp_a++;
        printf("\n");
    }
}
