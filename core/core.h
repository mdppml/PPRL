//
// Created by Mete Akgun on 28.12.21.
//

#ifndef CORE_H
#define CORE_H

#include "Party.h"
#include "../utils/test_functions.h"
#include <thread>
#include <mutex>
#include <bitset>

double recn_time = 0;
double mul_triple_gen = 0;
double mul_ef_calc = 0;
double mul_z_calc = 0;
double rec_calc = 0;
double rec_transfer = 0;
double rec_read = 0;

/**
 * Perform the truncation operation which we use to keep the number of fractional bit consistent after Multiply operation
 * @param proxy
 * @param z: value we want to Truncate
 * @return truncated z is returned
 */
uint64_t Truncate(Party *const proxy, uint64_t z, int shift = FRACTIONAL_BITS) {
    switch (proxy->GetPRole()) {
        case proxy1:
            z = ArithmeticShift(z, shift);
            break;
        case proxy2:
            z = -1 * ArithmeticShift(-1 * z, shift);
            break;
        case helper:
            break;
    }
    return z;
}

uint64_t *Reconstruct(Party *const proxy, const uint64_t *const a, size_t sz, uint64_t mask= RING_SIZE) {
    if (proxy->GetPRole() != helper) {
        uint64_t *b = new uint64_t[sz];
        if (proxy->GetPRole() == proxy1 ) {
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < sz; i++) {
                AddValueToCharArray(a[i], &ptr);
            }
            thread thr1 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer1(), sz * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
            thr1.join();
            thr2.join();

            ptr = proxy->GetBuffer2();
            for(size_t i = 0; i < sz; i++) {
                b[i] = ConvertToLong(&ptr);
            }

        } else if (proxy->GetPRole() == proxy2) {
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < sz; i++) {
                AddValueToCharArray(a[i], &ptr);
            }
            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer2(), sz * 8);
            thr1.join();
            thr2.join();
            ptr = proxy->GetBuffer2();
            for(size_t i = 0; i < sz; i++) {
                b[i] = ConvertToLong(&ptr);
            }
        }
        for(size_t i = 0; i < sz; i++) {
            b[i] = (a[i] + b[i]) & mask;
        }
        return b;
    } else {
        return nullptr;
    }
}

uint64_t Reconstruct(Party *const proxy, uint64_t a, uint64_t mask= RING_SIZE) {
    uint64_t *result_array = Reconstruct(proxy, &a, 1, mask);
    uint64_t result = result_array[0];
    delete[] result_array;
    return result;
}

uint64_t *ReconstructNarrow(Party *const proxy, const uint64_t *const a, size_t sz, uint64_t ringbits) {
    auto mask = (1<< ringbits)-1;
    auto bsz = (size_t)ceil(ringbits/8.0);
    uint64_t *b = new uint64_t[sz];
    if ( proxy->GetPRole() == proxy1 ) {
        unsigned char *ptr = proxy->GetBuffer1();
        WriteToBuffer(a, ptr, sz, bsz);
        thread thr1 = thread(Send,proxy->GetSocketP2(), proxy->GetBuffer1(), sz*bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP2(), proxy->GetBuffer2(), sz*bsz);
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();
        ReadBuffer(b, ptr, sz, bsz);
#pragma omp parallel for num_threads(4)
        for(size_t i = 0; i < sz; i++) {
            b[i] = (a[i] + b[i]) & mask;
        }
    } else if ( proxy->GetPRole() == proxy2) {
        unsigned char *ptr = proxy->GetBuffer1();
        WriteToBuffer(a, ptr, sz, bsz);
        thread thr1 = thread(Send,proxy->GetSocketP1(), proxy->GetBuffer1(), sz*bsz);
        thread thr2 = thread(Receive,proxy->GetSocketP1(), proxy->GetBuffer2(), sz*bsz);
        thr2.join();
        thr1.join();
        ptr = proxy->GetBuffer2();
        ReadBuffer(b, ptr, sz, bsz);
#pragma omp parallel for num_threads(4)
        for(size_t i = 0; i < sz; i++) {
            b[i] = (a[i] + b[i]) & mask;
        }
    }
    return b;
}

/**Reconstruct a secret shared 2D array.*/
uint64_t** Reconstruct(Party *const proxy, const uint64_t *const *const a, size_t n_row, size_t n_col) {
    if (proxy->GetPRole() != helper) {
        uint64_t **b = new uint64_t*[n_row];
        if (proxy->GetPRole() == proxy1) {
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < n_row; i++) {
                b[i] = new uint64_t[n_col];
                for( size_t j = 0; j < n_col; j++) {
                    AddValueToCharArray(a[i][j], &ptr);
                }
            }
            thread thr1 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer1(), n_row * n_col * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), n_row * n_col * 8);
            thr1.join();
            thr2.join();
            ptr = proxy->GetBuffer2();
            for(size_t i = 0; i < n_row; i++) {
                for(size_t j = 0; j < n_col; j++) {
                    b[i][j] = ConvertToLong(&ptr);
                }
            }
        } else if (proxy->GetPRole() == proxy2) {
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < n_row; i++) {
                for( size_t j = 0; j < n_col; j++) {
                    AddValueToCharArray(a[i][j], &ptr);
                }
            }
            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), n_row * n_col * 8);
            thread thr2 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer2(), n_row * n_col * 8);
            thr1.join();
            thr2.join();
            ptr = proxy->GetBuffer2();
            for(size_t i = 0; i < n_row; i++) {
                b[i] = new uint64_t[n_col];
                for(size_t j = 0; j < n_col; j++) {
                    b[i][j] = ConvertToLong(&ptr);
                }
            }
        }
        for(size_t i = 0; i < n_row; i++) {
            for(size_t j = 0; j < n_col; j++) {
                b[i][j] = (a[i][j] + b[i][j]);
            }
        }
        return b;
    } else {
        return nullptr;
    }
}

/**Reconstruct a secret shared 3D array.*/
uint64_t*** Reconstruct(Party *const proxy, const uint64_t *const *const *const a, size_t n_matrices, size_t n_row, size_t n_col) {
    if(proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t ***b = new uint64_t**[n_matrices];
        unsigned char *ptr = proxy->GetBuffer1();
        // add the values of a to the buffer and initialize b
        for(size_t n = 0; n < n_matrices; n++) {
            b[n] = new uint64_t * [n_row];
            for(size_t i = 0; i < n_row; i++) {
                b[n][i] = new uint64_t[n_col];
                for(size_t j = 0; j < n_col; j++) {
                    AddValueToCharArray(a[n][i][j], &ptr);
                }
            }
        }
        // determine the
        int *partner_socket;
        if (proxy->GetPRole() == proxy1) {
            partner_socket = proxy->GetSocketP2();
        }
        else {
            partner_socket = proxy->GetSocketP1();
        }
        thread thr1 = thread(Send, partner_socket, proxy->GetBuffer1(), n_matrices * n_row * n_col * 8);
        thread thr2 = thread(Receive, partner_socket, proxy->GetBuffer2(), n_matrices * n_row * n_col * 8);
        thr1.join();
        thr2.join();
        ptr = proxy->GetBuffer2();

        for(size_t n = 0; n < n_matrices; n++) {
            for(size_t i = 0; i < n_row; i++) {
                for(size_t j = 0; j < n_col; j++) {
                    b[n][i][j] = ConvertToLong(&ptr);
                }
            }
        }

        for(size_t n = 0; n < n_matrices; n++) {
            for(size_t i = 0; i < n_row; i++) {
                for(size_t j = 0; j < n_col; j++) {
                    b[n][i][j] += a[n][i][j];
                }
            }
        }

        return b;
    }
    return nullptr;
}


uint64_t Add(Party *const proxy, uint64_t a, uint64_t b) {
    return a + b;
}

/**
 * Adds values of a and b at equal position.
 * @param proxy
 * @param a
 * @param b
 * @param size length of vectors a and b
 * @return vector of length size containing the sum of according values in a and b.
 */
uint64_t* Add(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size) {
    uint64_t* sum = new uint64_t[size];
    for(size_t i = 0; i<size; i++){
        sum[i] = a[i] + b[i];
    }
    return sum;
}

/** For smaller ring size
 * @param mt1 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param mt2 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param size the number of multiplication triples that will be generated
 */
void GenerateMultiplicationTriple(Party *const proxy, uint64_t *const c1, size_t size, uint64_t mask) {

    for(size_t i = 0; i < size; i++) {
        uint64_t a0 = proxy->GenerateCommonRandom()&mask;
        uint64_t a1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t b0 = proxy->GenerateCommonRandom()&mask;
        uint64_t b1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t c0=  proxy->GenerateCommonRandom()&mask;
        c1[i] = (((a0+a1)*(b0+b1)) - c0)&mask;
        }
}

/**
 * Adds values of all vectors in a at equal position in a row to calculate their sum (sum over column where each row is one vector).
 * @param proxy
 * @param a matrix containing several vectors of length size. Values of all vectors at same position shall be summed up.
 * @param n_vectors number of vectors in a
 * @param size length of each vector in a
 * @return vector of length size
 */
uint64_t* Add(Party *const proxy, const uint64_t *const *const a, size_t n_vectors, size_t size) {
    uint64_t* res = new uint64_t [size];
    for(size_t i = 0; i<size; i++){
        res[i] = 0;
        for(size_t v = 0; v<n_vectors; v++){
            res[i] += a[v][i];
        }
    }
    return res;
}


/** Strictly a function for the helper to generate Beaver's multiplication triples.
 * It is used to compute proxy2's shares of c, which then need to be sent to it.
 * @param c1 array in which proxy2's shares of C will be stored
 * @param size the number of multiplication triples that will be generated
 */
 void GenerateMultiplicationTriple(Party *const proxy, uint64_t *const c1, size_t size) {
     for(size_t i = 0; i < size; i++) {
        uint64_t a0 = proxy->GenerateCommonRandom();
        uint64_t a1 = proxy->GenerateCommonRandom2();
        uint64_t b0 = proxy->GenerateCommonRandom();
        uint64_t b1 = proxy->GenerateCommonRandom2();
        uint64_t c0=  proxy->GenerateCommonRandom();
        c1[i] = ((a0+a1)*(b0+b1)) - c0;
    }
}

/** For smaller ring size and symmetric Multiply comm.
 * Half of the c0 and c1 is generated with RNG, rest is calculated and sent by Helper
 * @param mt1 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param mt2 3-by-size array whose rows will be a_i, b_i and c_i, respectively
 * @param size the number of multiplication triples that will be generated
 */
[[maybe_unused]] void GenerateMultiplicationTripleSym(Party *const proxy, uint64_t *const c0, uint64_t *const c1, size_t size, uint64_t mask) {
    for(size_t i = 0; i < size/2; i++) {
        uint64_t a0 = proxy->GenerateCommonRandom()&mask;
        uint64_t a1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t b0 = proxy->GenerateCommonRandom()&mask;
        uint64_t b1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t c0=  proxy->GenerateCommonRandom()&mask;
        c1[i] = (((a0+a1)*(b0+b1)) - c0)&mask; //(a0+a1)*(b0+b1) - c0
    }
    for(size_t i = 0; i < size/2; i++) {
        uint64_t a0 = proxy->GenerateCommonRandom()&mask;
        uint64_t a1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t b0 = proxy->GenerateCommonRandom()&mask;
        uint64_t b1 = proxy->GenerateCommonRandom2() & mask;
        uint64_t c1= proxy->GenerateCommonRandom2() & mask;
        c0[i] = (((a0+a1)*(b0+b1)) - c1)&mask; //(a0+a1)*(b0+b1) - c1
    }
}

/**
 * Adds values of all matrices in a at equal position to calculate their sum (sum over all matrices in a).
 * @param proxy
 * @param a 3-dmatrix containing several 2-d matrices in dimension rows x cols. Values of all matrices at same position shall be summed up.
 * @param n_matrices number of matrices in a
 * @param rows height of each matrix in a
 * @param cols width of each matrix in a
 * @return 2-d matrix of shape rows x cols with the summed up values.
 */
uint64_t** Add(Party *const proxy, const uint64_t *const *const *const a, size_t n_matrices, size_t rows, size_t cols) {
    uint64_t** res = new uint64_t *[rows];
    for(size_t r = 0; r<rows; r++){
        res[r] = new uint64_t [cols];
        for(size_t c = 0; c<cols; c++){
            res[r][c] = 0;
            for(size_t m = 0; m < n_matrices; ++m) {
                res[r][c] += a[m][r][c];
            }
        }
    }
    return res;
}
uint32_t* MultiplexNarrow(Party *const proxy, const uint32_t *const x, const uint32_t *const y, const uint32_t *const b, size_t sz, size_t ringbits) {
    auto mask = (1<< ringbits)-1;
    size_t bsz = ceil((double)ringbits/8.0);
    if (proxy->GetPRole() == proxy1){
        unsigned char *ptr = proxy->GetBuffer1();
        uint32_t *res = new uint32_t[sz];
        uint32_t *m1 = new uint32_t[sz];
        for (size_t i = 0; i < sz; i++) {
            uint32_t r1= proxy->GenerateCommonRandom();
            uint32_t r2= r1 & mask;
            r1 = (r1 >> ringbits) & mask;
            uint32_t r3= proxy->GenerateCommonRandom();
            uint32_t r4= r3 & mask;
            r3 = (r3 >> ringbits) & mask;

            m1[i] = ((b[i] * (x[i] - y[i])) - (r2*b[i]) - (r3*(x[i] - y[i])) - (r3*r4)) & mask ;
            uint32_t m2 = (b[i] + r1) & mask;
            uint32_t m3 = (x[i] - y[i] + r4) & mask;

            AddValueToCharArray(m2, &ptr, bsz);
            AddValueToCharArray(m3, &ptr, bsz);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * bsz * 2);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * bsz);
        ptr = proxy->GetBuffer1();
        for (size_t i = 0; i < sz; i++) {
            res[i] = (m1[i] + ConvertToLong(&ptr,bsz)) & mask;
            res[i] = (x[i] - res[i]) & mask;
        }
        delete [] m1;
        return res;

    }else if (proxy->GetPRole() == proxy2){
        unsigned char *ptr = proxy->GetBuffer1();
        uint32_t *res = new uint32_t[sz];
        uint32_t *m1 = new uint32_t[sz];
        for (size_t i = 0; i < sz; i++) {
            uint32_t r1= proxy->GenerateCommonRandom();
            uint32_t r2= r1 & mask;
            r1 = (r1 >> ringbits) & mask;
            uint32_t r3= proxy->GenerateCommonRandom();
            uint32_t r4= r3 & mask;
            r3 = (r3 >> ringbits) & mask;

            m1[i] = ((b[i] * (x[i] - y[i])) - (r1*(x[i] - y[i])) - (r1*r2) - (r4*b[i])) & mask;
            uint32_t m2 = (x[i] - y[i] + r2) & mask;
            uint32_t m3 = (b[i] + r3) & mask;

            AddValueToCharArray(m2, &ptr, bsz);
            AddValueToCharArray(m3, &ptr, bsz);

        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 2 * bsz);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * bsz);
        ptr = proxy->GetBuffer1();
        for (size_t i = 0; i < sz; i++) {
            res[i] = (m1[i] + ConvertToLong(&ptr,bsz)) & mask;
            res[i] = (x[i] - res[i]) & mask;
        }
        delete [] m1;
        return res;

    }else if (proxy->GetPRole() == helper){
        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 2 * bsz);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 2 * bsz);
        thr1.join();
        thr2.join();
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        uint32_t *m = new uint32_t[sz];
        for (size_t i = 0; i < sz; i++) {
            uint32_t m2 = ConvertToLong(&ptr, bsz);
            uint32_t m3 = ConvertToLong(&ptr, bsz);

            uint32_t m5 = ConvertToLong(&ptr2, bsz);
            uint32_t m6 = ConvertToLong(&ptr2, bsz);

            m[i] = ((m2 * m5) + (m3 * m6)) & mask;
        }
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (size_t i = 0; i < sz; i++) {
            uint32_t tmp1= proxy->GenerateRandom() & mask;
            uint32_t tmp2 = (m[i]-tmp1) & mask;
            AddValueToCharArray(tmp1, &ptr_out, bsz);
            AddValueToCharArray(tmp2, &ptr_out2, bsz);
        }
        delete []m;
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * bsz);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * bsz);
        thr1.join();
        thr2.join();
        return NULL;
    }
    return NULL;
}



uint64_t* Multiplex(Party *const proxy, const uint64_t *const x, const uint64_t *const y, const uint64_t *const b, size_t sz, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == proxy1){
        unsigned char *ptr = proxy->GetBuffer1();
        uint64_t *res = new uint64_t[sz];
        uint64_t *m1 = new uint64_t[sz];
        for (size_t i = 0; i < sz; i++) {
            uint64_t r1= proxy->GenerateCommonRandom(), r2= proxy->GenerateCommonRandom(), r3= proxy->GenerateCommonRandom(), r4= proxy->GenerateCommonRandom();

            m1[i] = (b[i] * (x[i] - y[i])) - (r2*b[i]) - (r3*(x[i] - y[i])) - (r3*r4);
            uint64_t m2 = b[i] + r1;
            uint64_t m3 = x[i] - y[i] + r4;

            AddValueToCharArray(m2, &ptr);
            AddValueToCharArray(m3, &ptr);
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        for (size_t i = 0; i < sz; i++) {
            res[i] = m1[i] + ConvertToLong(&ptr);
            res[i] = res[i] >> shift;
            res[i] = x[i] - res[i];
        }
        delete [] m1;
        return res;

    }else if (proxy->GetPRole() == proxy2){
        unsigned char *ptr = proxy->GetBuffer1();
        uint64_t *res = new uint64_t[sz];
        uint64_t *m1 = new uint64_t[sz];
        for (size_t i = 0; i < sz; i++) {
            uint64_t r1= proxy->GenerateCommonRandom(), r2= proxy->GenerateCommonRandom(), r3= proxy->GenerateCommonRandom(), r4= proxy->GenerateCommonRandom();

            m1[i] = (b[i] * (x[i] - y[i])) - (r1*(x[i] - y[i])) - (r1*r2) - (r4*b[i]);
            uint64_t m2 = x[i] - y[i] + r2;
            uint64_t m3 = b[i] + r3;

            AddValueToCharArray(m2, &ptr);
            AddValueToCharArray(m3, &ptr);

        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 16);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * 8);
        ptr = proxy->GetBuffer1();
        for (size_t i = 0; i < sz; i++) {
            res[i] = m1[i] + ConvertToLong(&ptr);
            res[i] = -1 * ((-1 * res[i]) >> shift);
            res[i] = x[i] - res[i];
        }
        delete [] m1;
        return res;

    }else if (proxy->GetPRole() == helper){
        thread thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 16);
        thread thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 16);
        thr1.join();
        thr2.join();
        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for (size_t i = 0; i < sz; i++) {
            uint64_t m2 = ConvertToLong(&ptr);
            uint64_t m3 = ConvertToLong(&ptr);

            uint64_t m5 = ConvertToLong(&ptr2);
            uint64_t m6 = ConvertToLong(&ptr2);

            uint64_t m = (m2 * m5) + (m3 * m6);
            m2 = proxy->GenerateRandom();
            m3 = m-m2;
            AddValueToCharArray(m2, &ptr_out);
            AddValueToCharArray(m3, &ptr_out2);
        }
        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 8);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 8);
        thr1.join();
        thr2.join();
        return NULL;
    }
    return NULL;
}

uint64_t Multiplex(Party *const proxy, uint64_t x, uint64_t y, uint64_t b, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Multiplex(proxy, nullptr, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t *result_vector = Multiplex(proxy, &x, &y, &b, 1, shift);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

/**
 * Private Compare Bool: Check b>a
 *
 * @param a reconstructed value
 * @param b boolean share
 * @param L1
 * @return
 */
uint8_t *PrivateCompareBool(Party *const proxy, const uint64_t *const a, const uint8_t *const b, size_t sz, int L1) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        for(size_t j = 0; j < sz; j++) {
            size_t jk = j * L1;
            uint8_t w_sum = 0;
            for (int i = L1 - 1; i >= 0; i--) {
                uint8_t a_bit = Bit(a[j], i);
                size_t k = jk + i;
                uint8_t w = Mod((b[k] + proxy->GetPRole() * a_bit - 2 * a_bit * b[k]) % LP, LP);
                proxy->GetBuffer1()[k] =
                        (Mod((proxy->GetPRole() * a_bit - b[k] + proxy->GetPRole() + w_sum), LP) * (proxy->GenerateCommonRandom() % (LP - 1) + 1)) %
                        LP;
                w_sum = (w_sum + w) % LP;
            }
            for (int i = 0; i < L1; i++) {
                size_t ind1 = (proxy->GenerateCommonRandom() % L1) + jk;
                size_t ind2 = (proxy->GenerateCommonRandom() % L1) + jk;
                uint8_t tmp = proxy->GetBuffer1()[ind1];
                proxy->GetBuffer1()[ind1] = proxy->GetBuffer1()[ind2];
                proxy->GetBuffer1()[ind2] = tmp;
            }
        }
        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * L1);
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz);
        uint8_t *r = new uint8_t[sz];
        for(size_t i = 0; i < sz; i++)
            r[i] = proxy->GetBuffer1()[i];
        return r;

    } else if (proxy->GetPRole() == helper) {
        Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), sz * L1);
        Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), sz * L1);
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for(size_t j = 0; j < sz; j++) {
            size_t jk = j * L1;
            uint8_t res = 0;
            for (int i = 0; i < L1; i++) {
                proxy->GetBuffer1()[jk + i] = (proxy->GetBuffer1()[jk + i] + proxy->GetBuffer2()[jk + i]) % LP;
                if (((int) proxy->GetBuffer1()[jk + i]) == 0) {
                    res = 1;
                    break;
                }
            }
            uint8_t res1 = proxy->GenerateRandom() % 2;
            uint8_t res2 = res ^res1;

            AddValueToCharArray(res1, &ptr_out);
            AddValueToCharArray(res2, &ptr_out2);
        }
        Send(proxy->GetSocketP1(), proxy->GetBuffer1(), sz);
        Send(proxy->GetSocketP2(), proxy->GetBuffer2(), sz);
        return nullptr;
    }
    return nullptr;
}

/** Check @p b>@p a
 *
 * @param a reconstructed value
 * @param b boolean share
 * @param L1
 * @return
 */
uint8_t PrivateCompareBool(Party *const proxy, uint64_t a, const uint8_t *const b, int L1) {
    if (proxy->GetPRole() == helper) {
        PrivateCompareBool(proxy, nullptr, nullptr, 1, L1);
        return 0;
    } else {
        auto result_vector = PrivateCompareBool(proxy, &a, b, 1, L1);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

/** Multiple modular conversions.
 *
 * @param x an array of values in the ring 2^63
 * @param sz the length of @p x
 * @return an array of values in the ring 2^64
 */
uint64_t *ModularConversion(Party *const proxy, const uint64_t *const x, size_t sz) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L_BIT - 1)];
        uint8_t *w = new uint8_t[sz];
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (8 + L_BIT));
        unsigned char *ptr = proxy->GetBuffer1();
        for(size_t i = 0; i < sz; i++) {
            ya[i] = ConvertToLong(&ptr);
            ConvertToArray(&ptr, &yb[i * (L_BIT - 1)], L_BIT - 1);
            w[i] = (*ptr);
            ptr++;
            z_1[i] = (x[i] + ya[i]) & N1_MASK;
        }
        uint64_t *z = Reconstruct(proxy, z_1, sz, N1_MASK);
        uint8_t *wc = PrivateCompareBool(proxy, z, yb, sz, L_BIT - 1);

        for(size_t i = 0; i < sz; i++) {
            w[i] = w[i] ^ wc[i];
            if (proxy->GetPRole() == proxy1 && z_1[i] > z[i])
                z_1[i] = z_1[i] + N1;
            z_1[i] = (z_1[i] - (ya[i] + w[i] * N1));
        }
        delete[] ya;
        delete[] yb;
        delete[] w;
        return z_1;
    }
    else if (proxy->GetPRole() == helper) {
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for(size_t i = 0; i < sz; i++) {
            uint64_t y = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            AddValueToCharArray(ya_1, &ptr_out);
            AddValueToCharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L_BIT - 1; j++) {
                uint8_t k = (y >> j) & 0x1;
                uint8_t yb_1 = proxy->GenerateRandom() % LP;
                uint8_t yb_2 = Mod(k - yb_1, LP);
                AddValueToCharArray(yb_1, &ptr_out);
                AddValueToCharArray(yb_2, &ptr_out2);
            }
            uint8_t w = 0;
            if (ya_1 > y)
                w = 1;
            uint8_t w_1 = proxy->GenerateRandom() % 2;
            uint8_t w_2 = w ^w_1;
            AddValueToCharArray(w_1, &ptr_out);
            AddValueToCharArray(w_2, &ptr_out2);
        }

        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (8 + L_BIT));
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (8 + L_BIT));
        thr1.join();
        thr2.join();
        // proxy1 and proxy2 will call MPrivateCompareBool
        PrivateCompareBool(proxy, 0, 0, sz, L_BIT - 1);
        return NULL;
    }
    return NULL;
}

/** Modular conversion.
 *
 * @param x a value in the ring 2^63
 * @return
 */
uint64_t ModularConversion(Party *const proxy, uint64_t x) {
    if (proxy->GetPRole() == helper) {
        ModularConversion(proxy, nullptr, 1);
        return 0;
    } else {
        uint64_t *result_vector = ModularConversion(proxy, &x, 1);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

void MostSignificantBitSubroutine(
    Party *const proxy,
    const uint64_t *const x,
    const uint64_t *const z,
    uint64_t *const z_1,
    const uint8_t *const yb,
    const uint64_t *const ya,
    uint8_t f,
    uint8_t rnd,
    size_t start_index,
    size_t end_index
){
    int L1 = L_BIT - 1;
    unsigned char *ptr_out = proxy->GetBuffer1();
    ptr_out += (start_index * (L1+16));
    size_t buffer_index = (start_index * (L1+16));
    size_t y_index = (start_index * L1);
    for(size_t i = start_index; i < end_index; i++) {
        uint8_t w_sum = 0;
        for (int t = L1 - 1; t >= 0; t--) {
            uint8_t a_bit = Bit(z[i], t);
            size_t bi = buffer_index + t;
            size_t yi = y_index + t;
            uint8_t w = Mod((yb[yi] + proxy->GetPRole() * a_bit - 2 * a_bit * yb[yi]) % LP, LP);
            proxy->GetBuffer1()[bi] = (Mod((proxy->GetPRole() * a_bit - yb[yi] + proxy->GetPRole() + w_sum), LP) * ((rnd % (LP - 1)) + 1)) % LP;
            rnd += 7;
            w_sum = (w_sum + w) % LP;
        }
        buffer_index += L1;
        y_index += L1;
        ptr_out += L1;


        uint8_t isWrap = 0;
        if (z[i]<z_1[i])
            isWrap = 1;
        z_1[i] =  z_1[i] + proxy->GetPRole() * isWrap * N1;
        AddValueToCharArray(proxy->GetPRole() * f * N1 - x[i] + z_1[i] - ya[i], &ptr_out);
        AddValueToCharArray(proxy->GetPRole() * (1 - f) * N1 - x[i] + z_1[i] - ya[i], &ptr_out);
        buffer_index +=16;
    }
}


/**
 * @brief Computes the most significant bit.
 *
 * MSB has 4 communication round. ModularConversion and PC are hardcoded in MostSignificantBit to reduce the number of
 * communication rounds of MostSignificantBit calls.
 * @param proxy p_proxy:
 * @param x p_x:
 * @param sz p_sz: The number of elements in x.
 * @param format p_format: Whether to convert the result to the regular representation. Defaults to true.
 * @return uint64_t* x < 0
 */
uint64_t *MostSignificantBit(Party *const proxy, const uint64_t *const x, size_t sz, int shift=FRACTIONAL_BITS, bool format = true) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint8_t f = proxy->GenerateCommonRandomByte() & 0x1;
        uint64_t *z_1 = new uint64_t[sz];
        uint64_t *ya = new uint64_t[sz];
        uint8_t *yb = new uint8_t[sz * (L_BIT - 1)];

        Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (8 + L_BIT - 1));

        unsigned char *ptr = proxy->GetBuffer1();
        for(size_t i = 0; i < sz; i++) {
            uint64_t dk = x[i] & N1_MASK;
            ya[i] = ConvertToLong(&ptr);
            ConvertToArray(&ptr, &yb[i * (L_BIT - 1)], L_BIT - 1);
            z_1[i] = (dk + ya[i]) & N1_MASK;
        }

        uint64_t *z = Reconstruct(proxy, z_1, sz, N1_MASK);
        size_t block_size = (size_t)ceil(sz * 1.0 / SOCKET_NUMBER);
        if (block_size == 0)
            block_size = sz;

        thread thr[SOCKET_NUMBER];
        size_t start_index = 0;
        size_t end_index = block_size;
        int thr_num = 0;
        for (int i = 0; i < SOCKET_NUMBER; i++) {
            uint8_t rnd = proxy->GenerateCommonRandomByte();
            thr[i] = thread(MostSignificantBitSubroutine, proxy, x, z, z_1, yb, ya, f, rnd, start_index, end_index);
            thr_num +=1;
            start_index += block_size;
            end_index += block_size;
            if (start_index >= sz)
                break;
            if (end_index > sz)
                end_index = sz;
        }
        for (int i = 0; i < thr_num; i++) {
            thr[i].join();
        }

        delete [] yb;
        delete [] ya;
        delete [] z_1;
        delete[] z;

        Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), sz * (16 + L_BIT - 1));
        Receive(proxy->GetSocketHelper(), proxy->GetBuffer2(), sz * 16);

        ptr = proxy->GetBuffer2();
        uint64_t *m = new uint64_t[sz];
        uint64_t val[2];
        for(size_t i = 0; i < sz; i++) {
            val[0] = ConvertToLong(&ptr);
            val[1] = ConvertToLong(&ptr);
            m[i] = val[f];
        }
        return m;

    }
    else if (proxy->GetPRole() == helper) {
        unsigned char *ptr_out = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        uint8_t *w = new uint8_t [sz];
        for(size_t i = 0; i < sz; i++) {
            uint64_t y = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_1 = proxy->GenerateRandom() & N1_MASK;
            uint64_t ya_2 = (y - ya_1) & N1_MASK;
            AddValueToCharArray(ya_1, &ptr_out);
            AddValueToCharArray(ya_2, &ptr_out2);
            for (int j = 0; j < L_BIT - 1; j++) {
                uint8_t k = (y >> j) & 0x1;
                uint8_t yb_1 = proxy->GenerateRandomByte() % 0x3f;
                uint8_t yb_2 = LP - yb_1 + k;
                AddValueToCharArray(yb_1, &ptr_out);
                AddValueToCharArray(yb_2, &ptr_out2);
            }
            w[i] = 0;
            if (y<ya_1)
                w[i] = 1;
        }
        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (8 + L_BIT - 1));
        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (8 + L_BIT - 1));
        thr1.join();
        thr2.join();

        thr1 = thread(Receive, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * (16 + L_BIT - 1));
        thr2 = thread(Receive, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * (16 + L_BIT - 1));
        thr1.join();
        thr2.join();


        unsigned char *ptr = proxy->GetBuffer1();
        unsigned char *ptr2 = proxy->GetBuffer2();
        ptr_out = proxy->GetBuffer1();
        ptr_out2 = proxy->GetBuffer2();

        int L1 = L_BIT-1;
        size_t jk = 0;
        for(size_t j = 0; j < sz; j++) {
            uint8_t res = 0;
            for (int i = 0; i < L1; i++) {
                uint8_t tmp = (proxy->GetBuffer1()[jk + i] + proxy->GetBuffer2()[jk + i]) % LP;
                if (((int) tmp) == 0) {
                    res = 1;
                    break;
                }
            }
            jk += L1;
            ptr += L1;
            ptr2 += L1;


            uint64_t val1 = (ConvertToLong(&ptr) + ConvertToLong(&ptr2)-(w[j]^res)*N1)/N1;
            uint64_t val2 = (ConvertToLong(&ptr) + ConvertToLong(&ptr2)-(w[j]^res)*N1)/N1;
            jk += 16;
            if(format) {
                val1 = ConvertToUint64((double) val1, shift);
                val2 = ConvertToUint64((double) val2, shift);
            }
            uint64_t vs_1 = proxy->GenerateRandom();
            uint64_t vs_2 = (val1 - vs_1);
            AddValueToCharArray(vs_1, &ptr_out);
            AddValueToCharArray(vs_2, &ptr_out2);
            vs_1 = proxy->GenerateRandom();
            vs_2 = (val2 - vs_1);
            AddValueToCharArray(vs_1, &ptr_out);
            AddValueToCharArray(vs_2, &ptr_out2);
        }

        thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), sz * 16);
        thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), sz * 16);
        thr1.join();
        thr2.join();

        delete [] w;

        return 0;
    }
    return 0;
}


/** Most significant bit: Returns the first (=left-most) bit of @p x.
 *
 * @param x
 * @return The first bit of @p x
 */
uint64_t MostSignificantBit(Party *const proxy, uint64_t x, int shift=FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        MostSignificantBit(proxy, nullptr, 1, shift);
        return 0;
    } else {
        auto result_array = MostSignificantBit(proxy, &x, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

uint64_t* ArithmeticAnd(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size, int shift=FRACTIONAL_BITS) {
    if (proxy->GetPRole() != helper) {
        std::unique_ptr<uint64_t[]> sums = std::make_unique<uint64_t[]>(size);
        uint64_t point_five = ConvertToUint64(0.5, shift);
        for(size_t i = 0; i < size; i++) {
            sums[i] = point_five - a[i] - b[i];
        }
        return MostSignificantBit(proxy, sums.get(), size, shift);
    } else {
        return MostSignificantBit(proxy, nullptr, size, shift);
    }
}

/**
 * @brief Compares all values of x to all values of y.
 *
 * In the resulting array, the first y_size results are all comparisons between the first x value and y and so on.
 * @param proxy p_proxy
 * @param x p_x: an array of shared values
 * @param y p_y:an array of shared values
 * @param x_size p_x_size: the amount of values in x
 * @param y_size p_y_size: the amount of values in y
 * @param shift p_shift: How many bits are reserved for the decimal places. Defaults to FRAC.
 * @return uint64_t* all comparison results between x and y
 */
uint64_t *CompareAll(
    Party *const proxy,
    const uint64_t *const x,
    const uint64_t *const y,
    size_t x_size,
    size_t y_size,
    int shift = FRACTIONAL_BITS
) {
    size_t overall_size = x_size*y_size;
    if (proxy->GetPRole() != helper) {
        uint64_t* difference = new uint64_t[overall_size];
        size_t overall_index;
        for(size_t x_index = 0; x_index < x_size; x_index++) {
            for(size_t y_index = 0; y_index < y_size; y_index++) {
                overall_index = x_index * y_size + y_index;
                difference[overall_index] = x[x_index] - y[y_index];
            }
        }
        uint64_t *result = MostSignificantBit(proxy, difference, overall_size, shift);
        delete[] difference;
        for(size_t i = 0; i < overall_size; i++) {
            result[i] = (proxy->GetPRole()<<shift) - (result[i] << shift);
        }
        return result;
    } else { // HELPER
        MostSignificantBit(proxy, nullptr, overall_size, shift);
        return nullptr;
    }
}

/** Comparison between two numbers.
 *
 * @param proxy
 * @param x
 * @param y
 * @return @p x > @p y
 */
uint64_t *Compare(Party *const proxy, const uint64_t *const x, const uint64_t *const y, size_t sz, int shift = FRACTIONAL_BITS) {
    if ( proxy->GetPRole() != helper) {
        uint64_t* diff = new uint64_t[sz];
        for (size_t i = 0; i < sz; i++) {
            diff[i] = y[i] - x[i];
        }
        uint64_t* m = MostSignificantBit(proxy, diff, sz, shift);
        delete[] diff;
        return m;
    }else if ( proxy->GetPRole() == helper) {
        MostSignificantBit(proxy, nullptr, sz, shift);
        return nullptr;
    }
    return nullptr;
}

/** Comparison between two numbers.
 *
 * @param proxy
 * @param x
 * @param y
 * @return @p x > @p y
 */
uint64_t Compare(Party *const proxy, uint64_t x, uint64_t y, int shift = FRACTIONAL_BITS) {
    if ( proxy->GetPRole() != helper) {
        uint64_t *result_array = Compare(proxy, &x, &y, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return  result;
    }else { // helper
        Compare(proxy, nullptr, nullptr, 1, shift);
        return 0;
    }
}


uint64_t* Equals(Party *const proxy, const uint64_t *const x, const uint64_t *const y, size_t size, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() != helper) {
        uint64_t *xyx = new uint64_t[size*3];
        std::copy(x, x+size, xyx);
        std::copy(y, y+size, xyx+size);
        std::copy(x, x+size, xyx+2*size);
        uint64_t* greater_and_smaller = Compare(proxy, xyx, xyx+size, size*2, shift);
        delete[] xyx;
        auto m = new uint64_t[size];
        for (size_t i = 0; i < size; i++) {
                m[i] = ConvertToUint64(0.5, shift) -greater_and_smaller[i] -greater_and_smaller[i+size];
        }
        delete[] greater_and_smaller;
        return m;
    } else {
        Compare(proxy, nullptr, nullptr, size*2, shift);
        return nullptr;
    }
}

uint64_t Equals(Party *const proxy, uint64_t x, uint64_t y, int shift = FRACTIONAL_BITS){
    if (proxy->GetPRole() == helper) {
        Equals(proxy, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t *result_array = Equals(proxy, &x, &y, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** Multiplication of two arrays of numbers.
 *
 * @param a one of the vectors of shares of the multiplicands
 * @param b the other vector of shares of the multiplicands
 * @param size the size of the vectors @p a and @p b
 * @return a vector containing the share of the result of the multiplication
 */
uint64_t *Multiply(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size, int shift = FRACTIONAL_BITS) {
    if (DEBUG_FLAG >= 1)
        cout << "************************************************************\nPMNF_MUL is called" << endl;
    if (proxy->GetPRole() == helper) {
        uint64_t *c1 = new uint64_t[size];

        GenerateMultiplicationTriple(proxy, c1, size);

        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for(size_t j = 0; j < size; j++) {
            AddValueToCharArray(c1[j], &ptr_out2);
        }

        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), size * 8);

        delete[] c1;

        if (DEBUG_FLAG >= 1)
            cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
        return 0;
    } else if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *mt[3];
        mt[0] = new uint64_t[size]; //a
        mt[1] = new uint64_t[size]; //b
        mt[2] = new uint64_t[size]; //c
        uint64_t *concat_e_f = new uint64_t[size * 2];
        if (proxy->GetPRole() == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * 8);
            unsigned char *ptr = proxy->GetBuffer1();
            for (size_t i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2();
                mt[1][i] = proxy->GenerateCommonRandom2();
                mt[2][i] = ConvertToLong(&ptr);
                concat_e_f[i] = a[i] - mt[0][i];
                concat_e_f[i + size] = b[i] - mt[1][i];
            }
        }
        else {
            for (size_t i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2();
                mt[1][i] = proxy->GenerateCommonRandom2();
                mt[2][i] = proxy->GenerateCommonRandom2();

                concat_e_f[i] = a[i] - mt[0][i];
                concat_e_f[i + size] = b[i] - mt[1][i];
            }
        }

        uint64_t *e_f = Reconstruct(proxy, concat_e_f, size * 2);
        uint64_t *e = e_f;
        uint64_t *f = &e_f[size];

        uint64_t *z = new uint64_t[size];
        for (size_t i = 0; i < size; i++) {
            z[i] = proxy->GetPRole() * e[i] * f[i] + f[i] * mt[0][i] + e[i] * mt[1][i] + mt[2][i];
            z[i] = Truncate(proxy, z[i], shift);
        }
        delete [] e_f;
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        if(DEBUG_FLAG >= 1)
            cout << "Returning from PMNF_MUL...\n************************************************************" << endl;
        return z;
    } else {
        return nullptr;
    }
}

 /** Multiplication of two numbers.
  *
  * @param proxy
  * @param a a share of the first multiplicand
  * @param b a share of the second multiplicand
  * @return the share of the multiplication of @p a and @p b
  */
uint64_t Multiply(Party *const proxy, uint64_t a, uint64_t b, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Multiply(proxy, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t *result_array = Multiply(proxy, &a, &b, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** Multiplication of two arrays of numbers for smaller ring size
 *
 * @param a one of the vectors of shares of the multiplicands
 * @param b the other vector of shares of the multiplicands
 * @param size the size of the vectors @p a and @p b
 * @return a vector containing the share of the result of the multiplication
 */
uint64_t *MultiplyNarrow(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size, size_t ringbits, int shift = FRACTIONAL_BITS) {
    uint64_t mask = (1<<ringbits)-1;
    size_t bsz = ceil((double)ringbits/8.0);
    if (proxy->GetPRole() == helper) {
        uint64_t *c1 = new uint64_t[size];

        GenerateMultiplicationTriple(proxy, c1, size, mask);
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        WriteToBuffer(c1, ptr_out2, size, bsz);
        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), size * bsz);

        delete[] c1;
        return 0;

    } else if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *mt[3];
        mt[0] = new uint64_t[size]; //a
        mt[1] = new uint64_t[size]; //b
        mt[2] = new uint64_t[size]; //c
        uint64_t *concat_e_f = new uint64_t[size * 2];

        if (proxy->GetPRole() == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), size * bsz);
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2() & mask;
                mt[1][i] = proxy->GenerateCommonRandom2() & mask;
                mt[2][i] = ConvertToLong(&ptr, bsz);
	    }

#pragma omp parallel for num_threads(4)
            for(size_t i = 0; i < size; ++i) {
                concat_e_f[i] = (a[i] - mt[0][i])&mask;
                concat_e_f[i + size] = (b[i] - mt[1][i])&mask;
            }
        }
        else { // P1
            for(size_t i = 0; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2() & mask;
                mt[1][i] = proxy->GenerateCommonRandom2() & mask;
                mt[2][i] = proxy->GenerateCommonRandom2() & mask;
		}

#pragma omp parallel for num_threads(4)
            for(size_t i = 0; i < size; ++i) {
                concat_e_f[i] = (a[i] - mt[0][i])&mask;
                concat_e_f[i + size] = (b[i] - mt[1][i])&mask;
             }
        }

        uint64_t *e_f = ReconstructNarrow(proxy, concat_e_f, size * 2, ringbits);
        uint64_t *e = e_f;
        uint64_t *f = &e_f[size];

        uint64_t *z = new uint64_t[size];
#pragma omp parallel for num_threads(4)
        for(size_t i = 0; i < size; i++) {
            z[i] = (proxy->GetPRole() * e[i] * f[i] + f[i] * mt[0][i] + e[i] * mt[1][i] + mt[2][i])&mask;
           // z[i] = Truncate(proxy, z[i], shift);
        }

        delete [] e_f;
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        return z;
    } else {
        return nullptr;
    }
}


/** Multiplication of two arrays of numbers for smaller ring size
 * Uses common randon generator in generation of beaver triples to reduce the communication cost
 * SYMMETRIC: Helper sending half of the c1 to P2  and half of the c0 to P1
 * @param a one of the vectors of shares of the multiplicands
 * @param b the other vector of shares of the multiplicands
 * @param size the size of the vectors @p a and @p b
 * @return a vector containing the share of the result of the multiplication
 */
 [[maybe_unused]] uint64_t *MUL_sym(
     Party *const proxy,
     const uint64_t *const a,
     const uint64_t *const b,
     size_t size,
     size_t ringbits,
     int shift = FRACTIONAL_BITS
 ) {
    uint64_t mask = (1<<ringbits)-1;
    size_t bsz = ceil((double)ringbits/8.0);
    size_t hsize = size/2; //half size
    if (proxy->GetPRole() == helper) {
        uint64_t *c0 = new uint64_t[hsize];
        uint64_t *c1 = new uint64_t[hsize];

        GenerateMultiplicationTripleSym(proxy,c0, c1, size, mask);
        unsigned char *ptr_out1 = proxy->GetBuffer1();
        unsigned char *ptr_out2 = proxy->GetBuffer2();
        for(size_t j = 0; j < hsize; j++) {
            AddValueToCharArray(c0[j], &ptr_out1, bsz);
            AddValueToCharArray(c1[j], &ptr_out2, bsz);
        }

        Send( proxy->GetSocketP1(), proxy->GetBuffer1(), hsize * bsz);
        Send( proxy->GetSocketP2(), proxy->GetBuffer2(), hsize * bsz);

        delete[] c0;
        delete[] c1;
        return 0;

    } else if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {

        uint64_t *mt[3];
        mt[0] = new uint64_t[size]; //a
        mt[1] = new uint64_t[size]; //b
        mt[2] = new uint64_t[size]; //c
        uint64_t *concat_e_f = new uint64_t[size * 2];

        if (proxy->GetPRole() == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), hsize * bsz);
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < hsize; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2() & mask;
                mt[1][i] = proxy->GenerateCommonRandom2() & mask;
                mt[2][i] = ConvertToLong(&ptr, bsz);

                concat_e_f[i] = (a[i] - mt[0][i])&mask;
                concat_e_f[i + size] = (b[i] - mt[1][i])&mask;
            }
            for(size_t i = hsize; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2() & mask;
                mt[1][i] = proxy->GenerateCommonRandom2() & mask;
                mt[2][i] = proxy->GenerateCommonRandom2() & mask;

                concat_e_f[i] = (a[i] - mt[0][i])&mask;
                concat_e_f[i + size] = (b[i] - mt[1][i])&mask;
            }
        }
        else { // P1
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), hsize * bsz);
            unsigned char *ptr = proxy->GetBuffer1();
            for(size_t i = 0; i < hsize; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2() & mask;
                mt[1][i] = proxy->GenerateCommonRandom2() & mask;
                mt[2][i] = proxy->GenerateCommonRandom2() & mask;

                concat_e_f[i] = (a[i] - mt[0][i])&mask;
                concat_e_f[i + size] = (b[i] - mt[1][i])&mask;
            }
            for(size_t i = hsize; i < size; ++i) {
                mt[0][i] = proxy->GenerateCommonRandom2() & mask;
                mt[1][i] = proxy->GenerateCommonRandom2() & mask;
                mt[2][i] = ConvertToLong(&ptr, bsz);

                concat_e_f[i] = (a[i] - mt[0][i])&mask;
                concat_e_f[i + size] = (b[i] - mt[1][i])&mask;
            }
        }
        uint64_t *e_f = ReconstructNarrow(proxy, concat_e_f, size * 2, ringbits);
        uint64_t *e = e_f;
        uint64_t *f = &e_f[size];

        uint64_t *z = new uint64_t[size];
        int pRole = proxy->GetPRole();
#pragma omp parallel for num_threads(4)
        for(size_t i = 0; i < size; i++) {
            z[i] = (pRole * e[i] * f[i] + f[i] * mt[0][i] + e[i] * mt[1][i] + mt[2][i])&mask;
            //z[i] = Truncate(proxy, z[i], shift);
        }
        delete [] e_f;
        delete [] concat_e_f;
        for (auto &i : mt) {
            delete[] i;
        }
        return z;
    } else {
        return nullptr;
    }
}

/** Multiple exponentials. Note that this function considers only the specific number of least significant bits not to
 * cause overflow. This is different for positive and negative powers.
 *
 * @param a the vector of values that will be used as the power of exp
 * @param size the length of @p a
 * @return a vector of arithmetic secret shares for each exp(@p a)
 */
uint64_t* Exp(Party *const proxy, const uint64_t *const a, size_t size, int shift = FRACTIONAL_BITS) {
    int p_role = proxy->GetPRole();
    int n_bits = proxy->GetNBits();
    int neg_n_bits = proxy->GetNegNBits();

    if (p_role == proxy1 || p_role == proxy2) {
        // compute the absolute of the input value
        uint64_t* msb_a = MostSignificantBit(proxy, a, size, shift);
        uint64_t* abs_a = new uint64_t[size];
        for(size_t i = 0; i < size; i++) {
            abs_a[i] = ((uint64_t) 0) - a[i];
        }

        // compute the possible contribution of positive and negative values
        uint64_t* pec = new uint64_t[n_bits];
        uint64_t* nec = new uint64_t[n_bits];
        if(p_role == proxy2) {
            for (int i = n_bits - 1; i >= 0; i--) {
                pec[n_bits - i - 1] = ConvertToUint64(exp(pow(2, i - shift)));
                if (i > neg_n_bits - 1) {
                    nec[n_bits - i - 1] = (((uint64_t) 1) << shift);
                } else {
                    nec[n_bits - i - 1] = ConvertToUint64(1.0 / exp(pow(2, i - shift)));
                }
            }
        }
        else {
            for (int i = n_bits - 1; i >= 0; i--) {
                pec[n_bits - i - 1] = 0;
                nec[n_bits - i - 1] = 0;
            }
        }

        // selection of the correct contribution from each bit of the input value based on the msb of the input value
        uint64_t *pos_e_contributions = new uint64_t[size * (n_bits + 1)]; // if the value is positive
        uint64_t *neg_e_contributions = new uint64_t[size * (n_bits + 1)]; // if the value is negative
        uint64_t *one_contributions = new uint64_t[size * n_bits]; // if the bit is zero regardless of the sign
        uint64_t *repeated_msb_a = new uint64_t[size * (n_bits + 1)]; // this will be used as a selection bit for the contributions of all bits

        for(size_t i = 0; i < size; i++) {
            pos_e_contributions[i * (n_bits + 1)] = a[i];
            neg_e_contributions[i * (n_bits + 1)] = abs_a[i];
            repeated_msb_a[i * (n_bits + 1)] = msb_a[i];
            for (int bi = 0; bi < n_bits; bi++) {
                pos_e_contributions[(i * (n_bits + 1)) + bi + 1] = pec[bi];
                neg_e_contributions[(i * (n_bits + 1)) + bi + 1] = nec[bi];
                one_contributions[(i * n_bits) + bi] = p_role * (((uint64_t) 1) << shift);
                repeated_msb_a[(i * (n_bits + 1)) + bi + 1] = msb_a[i];
            }
        }
        delete[] msb_a;
        delete[] abs_a;
        delete[] pec;
        delete[] nec;
        uint64_t *e_contributions = Multiplex(proxy, pos_e_contributions, neg_e_contributions, repeated_msb_a,
                                              size * (n_bits + 1), shift);
        delete[] pos_e_contributions;
        delete[] neg_e_contributions;
        delete[] repeated_msb_a;
        uint64_t* new_a = new uint64_t[size];
        uint64_t* selected_e_contributions = new uint64_t[size * n_bits];
        for(size_t i = 0; i < size; i++) {
            new_a[i] = e_contributions[i * (n_bits + 1)];
            for(size_t j = 0; j < n_bits; j++) {
                selected_e_contributions[(i * n_bits) + j] = e_contributions[(i * (n_bits + 1)) + j + 1];
            }
        }
        delete[] e_contributions;
        // arrange all the shifted versions of the input value for MostSignificantBit
        uint64_t *partial_a = new uint64_t[size * n_bits];
        for(size_t i = 0; i < size; i++) {
            for (size_t j = 0; j < n_bits; j++) {
                partial_a[(i * n_bits) + j] = new_a[i] << (L_BIT - n_bits + j);
            }
        }

        // get secret shared form of the bits of the values that could contribute into the result
        uint64_t *bit_shares = MostSignificantBit(proxy, partial_a, size * n_bits, shift);
        delete[] partial_a;
        // selection of the contribution of the bits of the value
        uint64_t *contributions = Multiplex(proxy, one_contributions, selected_e_contributions, bit_shares,
                                            size * n_bits, shift);
        delete[] one_contributions;
        delete[] bit_shares;
        // binary-tree-based multiplication of the contributions into the BenchmarkExp
        int cs = n_bits;
        bool flag = false;
        uint64_t* remaining = new uint64_t[size];
        uint64_t *tmp1, *tmp2;
        for(size_t j = 0; j < (size_t) ceil(log2(n_bits)); j++) {
            tmp1 = contributions;
            tmp2 = &contributions[cs / 2];

            if (cs % 2 == 1) {
                if (!flag) {
                    for(size_t i = 0; i < size; i++){
                        remaining[i] = contributions[(i * cs) + cs - 1];
                    }

                    tmp1 = new uint64_t[size * (cs / 2)];
                    tmp2 = new uint64_t[size * (cs / 2)];

                    for( size_t i = 0; i < size; i++) {
                        copy(contributions + (i * cs), contributions + (i * cs) + (cs / 2), tmp1 + (i * (cs / 2)));
                        copy(contributions + (i * cs) + (cs / 2), contributions + (i * cs) + 2 * (cs / 2), tmp2 + (i * (cs / 2)));
                    }

                    flag = true;
                } else {
                    tmp1 = new uint64_t[size * ((cs + 1) / 2)];
                    tmp2 = new uint64_t[size * ((cs + 1) / 2)];

                    for(size_t i = 0; i < size; i++) {
                        copy(contributions + (i * cs), contributions + (i * cs) + ((cs + 1) / 2), tmp1 + (i * ((cs + 1) / 2)));
                        copy(contributions + (i * cs) + ((cs + 1) / 2), contributions + ((i + 1) * cs), tmp2 + (i * ((cs + 1) / 2)));
                        tmp2[(i + 1) * ((cs + 1) / 2) - 1] = remaining[i];
                    }

                    cs++;
                    flag = false;
                }
            }
            else {
                tmp1 = new uint64_t[size * (cs / 2)];
                tmp2 = new uint64_t[size * (cs / 2)];

                for( size_t i = 0; i < size; i++) {
                    copy(contributions + (i * cs), contributions + (i * cs) + (cs / 2), tmp1 + (i * (cs / 2)));
                    copy(contributions + (i * cs) + (cs / 2), contributions + ((i + 1) * cs), tmp2 + (i * (cs / 2)));
                }
            }
            delete [] contributions;
            contributions = Multiply(proxy, tmp1, tmp2, size * (cs / 2), shift);

            delete [] tmp1;
            delete [] tmp2;

            cs /= 2;
        }

        // deleting dynamically allocated arrays
        delete [] remaining;

        return contributions;
    }
    else if ( p_role == helper) {
        MostSignificantBit(proxy, 0, size, shift);
        Multiplex(proxy, 0, 0, 0, size * (n_bits + 1), shift);
        MostSignificantBit(proxy, 0, size * n_bits, shift);
        Multiplex(proxy, 0, 0, 0, size * n_bits, shift);
        size_t current_size = n_bits;
        bool flag = false;
        for(size_t i = 0; i < (size_t) ceil(log2(n_bits)); i++) {
            if (current_size % 2 == 1) {
                if (!flag) {
                    flag = true;
                } else {
                    current_size++;
                    flag = false;
                }
            }
            Multiply(proxy, 0, 0, size * (current_size / 2), shift);
            current_size /= 2;
        }

        return 0;
    }
    else {
        return nullptr;
    }
}

/** Exponential. Note that this function considers only the specific number of least significant bits not to cause
 * overflow. This is different for positive and negative powers.
 *
 * @param a the value that will be used as the power of exp
 * @return Returns the arithmetic secret share of exp(@p a)
 */
uint64_t Exp(Party *const proxy, uint64_t a, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Exp(proxy, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t* result_array = Exp(proxy, &a, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** PartialSum: sum the elements of each section separately.
 *
 * @param a the vector on which we perform the partial summation
 * @param size the size of @p a
 * @param d the size of the part that we will use to Partition @p a
 * @return
 */
uint64_t* PartialSum(Party *const proxy, const uint64_t *const a, size_t size, uint32_t d) {
    int p_role = proxy->GetPRole();
    if(p_role == proxy1 || p_role == proxy2) {
        uint64_t *ps_x = new uint64_t[size / d];
        for(uint32_t base = 0; base < size; base += d) {
            uint64_t tmp = 0;
            for(uint32_t i = 0; i < d; i++) {
                tmp += a[base + i];
            }
            ps_x[base / d] = tmp;
        }
        return ps_x;
    }
    else {
        return NULL;
    }
}

/** Computes the dot product of arithmetically shared vectors, which are formed by vectors.
 *
 * @param a vector formed by vectors of given size
 * @param b vector formed by vectors of given size
 * @param size the length of the vectors
 * @param d the size of the partial vectors forming the main vectors
 * @return Dot product of the given (@p size / @p d) vectors as a vector of (@p size / @p d)
 */
uint64_t* DotProduct(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size, uint32_t d, int shift = FRACTIONAL_BITS) {
    int p_role = proxy->GetPRole();
    if(p_role == proxy1 || p_role == proxy2) {
        // compute elementwise multiplication of vectors
        uint64_t *ew_xy = Multiply(proxy, a, b, size, shift);
        // sum the vectors in the main vector
        uint64_t *dp_shr = PartialSum(proxy, ew_xy, size, d);

        delete [] ew_xy;

        return dp_shr;
    }
    else if(p_role == helper) {
        Multiply(proxy, 0, 0, size, shift);
        return NULL;
    }
    else {
        return NULL;
    }

}

/** computes the dot product of two single arithmetically shared vectors.
 *
 * @param proxy
 * @param a vector
 * @param b vector
 * @param size the length of the vectors
 * @return
 */
uint64_t DotProduct(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        DotProduct(proxy, nullptr, nullptr, size, size, shift);
        return 0;
    } else {
        uint64_t *result_vector = DotProduct(proxy, a, b, size, size, shift);
        uint64_t result = result_vector[0];
        delete[] result_vector;
        return result;
    }
}

/** Perform several multiplications of matrices of size a_row-by-a_col and a_col-by-b_col stored in a and b.
 *
 * @param a three dimensional matrix
 * @param b three dimensional matrix
 * @param n_matrices number of two-dimensional matrices of @p a and @p b
 * @param a_row number of rows per two-dimensional matrix
 * @param a_col number of columns per two-dimensional matrix of @p a
 * @param b_col number of columns per two-dimensional matrix of @p b
 * @return a matrix of size @p n_matrices by @p a_row by @p b_col
 */
uint64_t*** MatrixMatrixMultiply(
    Party *const proxy,
    const uint64_t *const *const *const a,
    const uint64_t *const *const *const b,
    size_t n_matrices,
    size_t a_row,
    size_t a_col,
    size_t b_col,
    int shift= FRACTIONAL_BITS
) {
    int p_role = proxy->GetPRole();
    if (p_role == proxy1 || p_role == proxy2) {
        // templates for E and F matrices
        uint64_t ***E = new uint64_t**[n_matrices];
        uint64_t ***F = new uint64_t**[n_matrices];

        // receive the shares of A, B and C matrices
        if (p_role == proxy2) {
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(),
                    n_matrices * a_row * b_col * 8);
        }
        unsigned char *ptr = proxy->GetBuffer1();

        uint64_t ***mt2 = new uint64_t**[n_matrices];

        for(size_t i = 0; i < n_matrices; i++) {
            // <X-A>_1
            E[i] = new uint64_t *[a_row];
            for(size_t j = 0; j < a_row; j++) {
                E[i][j] = new uint64_t[a_col];
                for(size_t k = 0; k < a_col; k++) {
                    E[i][j][k] = a[i][j][k] - proxy->GenerateCommonRandom2();
                }
            }

            // <Y-B>_1
            F[i] = new uint64_t *[a_col];
            for(size_t j = 0; j < a_col; j++) {
                F[i][j] = new uint64_t[b_col];
                for(size_t k = 0; k < b_col; k++) {
                    F[i][j][k] = b[i][j][k] - proxy->GenerateCommonRandom2();
                }
            }

            if (p_role == proxy1) { // proxy1 generates its share of C locally using the common random generator with the helper
                mt2[i] = new uint64_t *[a_row];
                for(size_t j = 0; j < a_row; j++) {
                    mt2[i][j] = new uint64_t[b_col];
                    for(size_t k = 0; k < b_col; k++) {
                        mt2[i][j][k] = proxy->GenerateCommonRandom2();
                    }
                }
            }
            else { // proxy2 gets its share of C from the helper
                ConvertTo2dArray(&ptr, mt2[i], a_row, b_col);
            }
        }

        // reconstruct E = X-A and F = Y-B
        uint64_t ***recE = Reconstruct(proxy, E, n_matrices, a_row, a_col);
        uint64_t ***recF = Reconstruct(proxy, F, n_matrices, a_col, b_col);

        for(size_t n = 0; n < n_matrices; n++) {
            for(size_t r = 0; r < a_row; r++) {
                delete [] E[n][r];
            }
            delete [] E[n];

            for(size_t r = 0; r < a_col; r++) {
                delete [] F[n][r];
            }
            delete [] F[n];
        }
        delete [] E;
        delete [] F;

        // compute -role * E * F + <X>_role * F + E * <Y>_role + <C>_role
        uint64_t ***z = new uint64_t**[n_matrices];
        uint64_t ***tmpEF;
        uint64_t ***tmpXF;
        uint64_t ***tmpEY;
        // the part below takes the majority of the time - it can be addressed for optimization
        tmpEF = LocalMatrixMatrixMultiply(recE, recF, n_matrices, a_row, a_col, b_col, 0);
        tmpXF = LocalMatrixMatrixMultiply(a, recF, n_matrices, a_row, a_col, b_col, 0);
        tmpEY = LocalMatrixMatrixMultiply(recE, b, n_matrices, a_row, a_col, b_col, 0);

        for(size_t n = 0; n < n_matrices; n++) {
            z[n] = new uint64_t * [a_row];
            for(size_t i = 0; i < a_row; i++) {
                z[n][i] = new uint64_t[b_col];
                for(size_t j = 0; j < b_col; j++) {
                    z[n][i][j] = (p_role * -1) * tmpEF[n][i][j] + tmpXF[n][i][j] + tmpEY[n][i][j] + mt2[n][i][j];
                    z[n][i][j] = Truncate(proxy, z[n][i][j], shift);
                }
            }
        }

        for(size_t n = 0; n < n_matrices; n++) {
            for(size_t r = 0; r < a_row; r++) {
                delete [] recE[n][r];
                delete [] tmpEF[n][r];
                delete [] tmpXF[n][r];
                delete [] tmpEY[n][r];
                delete [] mt2[n][r];
            }
            delete [] recE[n];
            delete [] tmpEF[n];
            delete [] tmpXF[n];
            delete [] tmpEY[n];
            delete [] mt2[n];

            for(size_t r = 0; r < a_col; r++) {
                delete [] recF[n][r];
            }
            delete [] recF[n];
        }
        delete [] recE;
        delete [] recF;
        delete [] tmpEF;
        delete [] tmpXF;
        delete [] tmpEY;
        delete [] mt2;

        return z;
    }
    else if(p_role == helper) {
        unsigned char *ptr_out = proxy->GetBuffer1();

        // temporary matrices to hold the current A, B and C matrices
        uint64_t **tmpA = new uint64_t*[a_row];
        uint64_t **tmpA1 = new uint64_t*[a_row];
        uint64_t **tmpA2 = new uint64_t*[a_row];
        uint64_t **tmpB = new uint64_t*[a_col];
        uint64_t **tmpB1 = new uint64_t*[a_col];
        uint64_t **tmpB2 = new uint64_t*[a_col];
        uint64_t **tmpC1 = new uint64_t*[a_row];
        uint64_t **tmpC;

        for(size_t i = 0; i < a_row; i++) {
            tmpA[i] = new uint64_t[a_col];
            tmpA1[i] = new uint64_t[a_col];
            tmpA2[i] = new uint64_t[a_col];
            tmpC1[i] = new uint64_t[b_col];
        }

        for(size_t i = 0; i < a_col; i++) {
            tmpB[i] = new uint64_t[b_col];
            tmpB1[i] = new uint64_t[b_col];
            tmpB2[i] = new uint64_t[b_col];
        }

        uint64_t tmp; // to hold the generated random values

        // matrix generations
        for(size_t n = 0; n < n_matrices; n++) {
            // generation of A and its shares
            for(size_t i = 0; i < a_row; i++) {
                for(size_t j = 0; j < a_col; j++) {
                    tmpA1[i][j] = proxy->GenerateCommonRandom();
                    tmpA2[i][j] = proxy->GenerateCommonRandom2();
                    tmpA[i][j] = tmpA1[i][j] + tmpA2[i][j];
                }
            }

            // generation of B and its shares
            for(size_t i = 0; i < a_col; i++) {
                for(size_t j = 0; j < b_col; j++) {
                    tmpB1[i][j] = proxy->GenerateCommonRandom();
                    tmpB2[i][j] = proxy->GenerateCommonRandom2();
                    tmpB[i][j] = tmpB1[i][j] + tmpB2[i][j];
                }
            }

            // calculation of A * B
            tmpC = LocalMatrixMatrixMultiply(tmpA, tmpB, a_row, a_col, b_col, 0); // why shift=0?

            // generation of shares of C
            for(size_t i = 0; i < a_row; i++) {
                for(size_t j = 0; j < b_col; j++) {
                    tmp = proxy->GenerateCommonRandom();
                    AddValueToCharArray(tmpC[i][j] - tmp, &ptr_out);
                }
            }

            for(size_t i = 0; i < a_row; i++) {
                delete [] tmpC[i];
            }
            delete [] tmpC;
        }

        for(size_t i = 0; i < a_row; i++) {
            delete [] tmpA[i];
            delete [] tmpA1[i];
            delete [] tmpA2[i];
            delete [] tmpC1[i];
        }
        delete [] tmpA;
        delete [] tmpA1;
        delete [] tmpA2;
        delete [] tmpC1;
        for(size_t i = 0; i < a_col; i++) {
            delete [] tmpB[i];
            delete [] tmpB1[i];
            delete [] tmpB2[i];
        }
        delete [] tmpB;
        delete [] tmpB1;
        delete [] tmpB2;

        // send these matrices to proxy1 and proxy2, respectively
        Send(proxy->GetSocketP2(), proxy->GetBuffer1(), n_matrices * a_row * b_col * 8);

        return nullptr;
    }
    else {
        return nullptr;
    }
}

/** Perform multiplication of matrices a and b.
 * The function assumes that the number of columns of a equals to the number of rows of b.
 *
 * @param a two dimensional matrix
 * @param b two dimensional matrix
 * @param a_row number of rows of @p a and @p b
 * @param a_col number of columns of @p a
 * @param b_col number of columns of @p b
 * @return a matrix of size @p a_row by @p b_col
 */
uint64_t** MatrixMatrixMultiply(
    Party *const proxy,
    const uint64_t *const *const a,
    const uint64_t *const *const b,
    size_t a_row,
    size_t a_col,
    size_t b_col,
    int shift = FRACTIONAL_BITS
) {
    if (proxy->GetPRole() == helper) {
        MatrixMatrixMultiply(proxy, nullptr, nullptr, 1, a_row, a_col, b_col, shift);
        return nullptr;
    } else {
        uint64_t ***result_array = MatrixMatrixMultiply(proxy, &a, &b, 1, a_row, a_col, b_col, shift);
        uint64_t **result = result_array[0];
        delete[] result_array;
        return result;
    }
}

/** Perform n_matrices multiplications of matrices of size a_row-by-a_col and vectors of size a_col stored in a and
 * b, respectively.
 *
 * @param a three dimensional matrix
 * @param b two dimensional matrix
 * @param n_matrices number of matrices in @p a / vectors in @p b
 * @param a_row number of rows of @p a
 * @param a_col number of columns of @p a / size of @p b
 * @return a two-dimensional matrix of size @p n_matrices by @p a_row
 */
uint64_t** MatrixVectorMultiply(
    Party *const proxy,
    const uint64_t *const *const *const a,
    const uint64_t *const *const b,
    size_t n_matrices,
    size_t a_row,
    size_t a_col,
    int shift = FRACTIONAL_BITS
) {
    int p_role = proxy->GetPRole();
    if (p_role == proxy1 || p_role == proxy2) {
        // reformat the given vectors into matrices, ensuring that one of their dimensions is 1
        uint64_t ***mat_b = new uint64_t**[n_matrices];
        for(size_t i = 0; i < n_matrices; i++) {
            mat_b[i] = new uint64_t *[a_col];
            for(size_t j = 0; j < a_col; j++) {
                mat_b[i][j] = new uint64_t[1];
                mat_b[i][j][0] = b[i][j];
            }
        }

        uint64_t ***result_array = MatrixMatrixMultiply(proxy, a, mat_b, n_matrices, a_row, a_col, 1);

        for(size_t i = 0; i < n_matrices; i++) {
            for(size_t j = 0; j < a_col; j++) {
                delete [] mat_b[i][j];
            }
            delete [] mat_b[i];
        }
        delete [] mat_b;

        uint64_t **result = new uint64_t *[n_matrices];
        for(uint64_t i = 0; i < n_matrices; i++) {
            result[i] = new uint64_t[a_row];
            for(uint64_t j = 0; j < a_row; j++) {
                result[i][j] = result_array[i][j][0];
            }
        }
        return result;
    }
    else if(p_role == helper) {
        MatrixMatrixMultiply(proxy, NULL, NULL, n_matrices, a_row, a_col, 1, shift);
        return NULL;
    }
    else {
        return nullptr;
    }
}

/** Perform multiplication of matrix a and vector b. The function assumes that the number of columns of a is equal to
 * the length of b.
 * @param a two dimensional matrix
 * @param b vector
 * @param a_row number of rows of @p a
 * @param a_col number of columns of @p a / size of @p b
 * @return a vector of size @p a_row
 */
uint64_t* MatrixVectorMultiply(
    Party *const proxy,
    const uint64_t *const *const a,
    const uint64_t *const b,
    size_t a_row,
    size_t a_col,
    int shift = FRACTIONAL_BITS
) {
    if (proxy->GetPRole() == helper) {
        MatrixVectorMultiply(proxy, nullptr, nullptr, 1, a_row, a_col, shift);
        return nullptr;
    } else {
        uint64_t** result_array = MatrixVectorMultiply(proxy, &a, &b, 1, a_row, a_col, shift);
        uint64_t* result = result_array[0];
        delete[] result_array;
        return result;
    }
}

// TODO test
/** Get the Modular Inverse (ModularInverse) of a given number a with modulo being the specified ring size.
 * For the resulting/returned value b, it must hold ab Mod(modulo) are congruent to 1. The modulo under which a
 * multiplied with the inverse are equal to 1, will always be the ring size.
 * @param a secret share of the value for which the modular inverse shall be calculated.
 * @return the secret share of the modular inverse of a under the ring size.
 */
uint64_t ModularInverse(Party *const proxy, uint64_t a){
    cout << "searching for ModularInverse of value " << ConvertToDouble(a) << endl;
    uint64_t exchangingBit = RING_SIZE / 64;
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        for (uint16_t step = 0; step < 64; step++) {
            cout << "step " << step << endl;
            uint64_t ringProducts [exchangingBit];
            // start with 1 because 0 does not have an inverse value.
            for (uint64_t x = 1; x <= exchangingBit; x++) {
                uint64_t modInv = x*step + x;
                ringProducts[x - 1] = (a * modInv) & RING_SIZE; // ModularConversion(proxy, t); ?
            }
            cout << "stored all ring products..." << endl;
            unsigned char *ptr = proxy->GetBuffer1();
            AddValueToCharArray(ringProducts, &ptr, exchangingBit);
            Send(proxy->GetSocketHelper(), proxy->GetBuffer1(), exchangingBit * 8);

            cout << "sent ring products to helper" << endl;
            // receive fresh share from helper
            Receive(proxy->GetSocketHelper(), proxy->GetBuffer1(), 8); // either the share of modInv or -1 to identify, to continue searching
            ptr = proxy->GetBuffer1();
            uint64_t share = ConvertToLong(&ptr);
            cout << "got fresh share from helper: " << share << endl;
            if (share != -1){
                // the modInv has been found
                cout << "ModularInverse was found: " << share << endl;
                return share;
            }
        }
        return 0;
    }
    else if (proxy->GetPRole() == helper) {
        for (uint16_t step = 0; step < 64; step++) {
            cout << "step " << step << endl;
            Receive(proxy->GetSocketP1(), proxy->GetBuffer1(), exchangingBit * 8);
            Receive(proxy->GetSocketP2(), proxy->GetBuffer2(), exchangingBit * 8);
            unsigned char *ptr1 = proxy->GetBuffer1();
            unsigned char *ptr2 = proxy->GetBuffer2();
            cout << "got ring products from parties..." << endl;

            uint64_t ringProducts_recon[exchangingBit];
            ringProducts_recon[0] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2)); //modInv = exchangeBit*step + 1
            uint64_t m;
            for(uint64_t i = 1; i < exchangingBit; i++){
                // reconstructed product of a was: exchangeBit * step + i+1
                ringProducts_recon[i] = (ConvertToLong(&ptr1) + ConvertToLong(&ptr2));
                for(uint64_t j = 0; j < i; j++){
                    if(((ringProducts_recon[j] + ringProducts_recon[i]) & RING_SIZE) == 1){
                        //Mod inverse of a is found: i+1 + j+1
                        m = exchangingBit * 2 * step + i + j + 2;
                        cout << "ModularInverse was found: " << m << endl;
                        // SEND fresh share of found modular inverse
                        //reassign buffer because ptr1 and ptr2 were incremented by convert2Long calls.
                        ptr1 = proxy->GetBuffer1();
                        ptr2 = proxy->GetBuffer2();

                        uint64_t tmp = proxy->GenerateRandom();
                        AddValueToCharArray(tmp, &ptr1);
                        AddValueToCharArray(m - tmp, &ptr2);

                        thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 8);
                        thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 8);
                        thr1.join();
                        thr2.join();
                        cout << "sent fresh share of ModularInverse to parties; m= " << m << endl;
                        return 0;
                    }
                }
            }
            //reassign buffer because ptr1 and ptr2 were incremented by convert2Long calls.
            ptr1 = proxy->GetBuffer1();
            ptr2 = proxy->GetBuffer2();

            uint64_t noValFound = -1;
            AddValueToCharArray(noValFound, &ptr1);
            AddValueToCharArray(noValFound, &ptr2);

            thread thr1 = thread(Send, proxy->GetSocketP1(), proxy->GetBuffer1(), 8);
            thread thr2 = thread(Send, proxy->GetSocketP2(), proxy->GetBuffer2(), 8);
            thr1.join();
            thr2.join();
        }
        return 0;
    }
    return -1;
}


// TODO test
/** Compute the vectorized division of a / b where a and b are vectors - not the integer approximation of the result
 *
 * @param proxy : Party instance
 * @param a : vector of dividends
 * @param b : vector of dividers
 * @param size: number of division operations - which is the size of a and b
 * @param first_call : indicates whether the aucDivide call is for the integer part of the division result, i.e. first call.
 * If it is the first call, then there will be the second call of aucDivide for the fractional part of the division result
 * @return vector (a / b)
 */
uint64_t* Divide(Party *const proxy, const uint64_t *a, const uint64_t *b, size_t size, int shift = FRACTIONAL_BITS, bool first_call = true) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *signs;
        if(first_call) {
            uint64_t *inp1 = new uint64_t[2 * size];
            uint64_t *inp2 = new uint64_t[2 * size];
            for(size_t i = 0; i < size; i++) {
                inp1[i] = a[i];
                inp1[size + i] = b[i];
                inp2[i] = (uint64_t) 0 - a[i];
                inp2[size + i] = (uint64_t) 0 - b[i];
            }
            signs = MostSignificantBit(proxy, inp1, 2 * size, shift);
            uint64_t *abs_vals = Multiplex(proxy, inp1, inp2, signs, 2 * size, shift);
            a = &abs_vals[0];
            b = &abs_vals[size];

            delete [] inp1;
            delete [] inp2;
        }

        // initialize the variables for quotient and remainder and the Role vector, which is a vector full of the Role value
        uint64_t *Q = new uint64_t[size];
        uint64_t *R = new uint64_t[size];

        // obtain every bit of the dividend
        uint64_t *msb_bits_of_a = new uint64_t[L_BIT * size];

        for(size_t i = 0; i < size; i++) { // each value
            Q[i] = 0;
            R[i] = 0;
            uint64_t tmp = a[i];
            for(int j = 0; j < L_BIT; j++) { // each bit of the value
                msb_bits_of_a[i * L_BIT + j] = tmp;
                tmp = tmp << 1;
            }
        }
        uint64_t *bits_of_a = MostSignificantBit(proxy, msb_bits_of_a, L_BIT * size, shift, false);

        delete [] msb_bits_of_a;

        // traverse all bits of the dividend
        for (int16_t j = L_BIT - 1; j >= 0; j--) {
            for(size_t i = 0; i < size; i++) {
                R[i] = R[i] << 1; // shift the remainder
                R[i] += bits_of_a[(i * L_BIT) + (L_BIT - 1 - j)];
            }

            uint64_t *c = Compare(proxy, R, b, size, shift); // compare the current R and divider

            uint64_t *o1 = new uint64_t[2 * size];
            uint64_t *o2 = new uint64_t[2 * size];
            for(size_t i = 0; i < size; i++) {
                o1[2 * i] = c[i];
                o1[2 * i + 1] = c[i];
                o2[2 * i] = b[i];
                o2[2 * i + 1] = ((uint64_t) proxy->GetPRole()) << j;
            }

            // if the current R is larger than or equal to the divider, subtract the divider from R
            uint64_t *v = Multiply(proxy, o1, o2, 2 * size, shift);
            for(size_t i = 0; i < size; i++) {
                R[i] = R[i] - v[2 * i];
                Q[i] = Q[i] + v[2 * i + 1];
            }
            delete [] c;
            delete [] o1;
            delete [] o2;
            delete [] v;
        }

        delete [] bits_of_a;

        if(first_call) {
            // determine the selection bits for the signs of the results based on the signs of a's and b's
            // choose the positive result if a < 0 and b < 0, or a >= 0 and b >= 0
            // choose the negative result if a >= 0 and b < 0, or a < 0 and b >= 0
            // This is exactly what XOR does. We mimic XOR arithmetically, i.e. a XOR b = a + b - 2ab
            uint64_t *tmp = Multiply(proxy, signs, &signs[size], size, shift); // for determining the signs of the results
            uint64_t *c = new uint64_t[size]; // for determining the signs of the results - selection bits
            for(size_t i = 0; i < size; i++) {
                R[i] = R[i] << shift; // prepare the remainder for the second division call
                Q[i] = Q[i] << shift; // prepare the quotient for the final quotient
                c[i] = (signs[i] + signs[i + size]) - 2 * tmp[i]; // for determining the signs of the results
            }
            delete [] tmp;

            uint64_t *neg_Q = new uint64_t[size]; // the negative of the results in case they are the correct ones based on signs
            uint64_t *sec_div = Divide(proxy, R, b, size, shift, false); // second division call for the fractional part of the final quotient
            for(size_t i = 0; i < size; i++) {
                Q[i] += sec_div[i];
                neg_Q[i] = (uint64_t) 0 - Q[i];
            }

            delete [] sec_div;
            delete [] signs;

            // based on the above analysis, we select the correct version of the final quotient
            uint64_t *div_res = Multiplex(proxy, Q, neg_Q, c, size, shift);
            delete [] c;
            delete [] neg_Q;
            delete [] Q;
            delete [] R;
            delete [] a;
            return div_res;
        }

        delete [] R;

        return Q;
    }
    else if (proxy->GetPRole() == helper) {
        if(first_call) {

            MostSignificantBit(proxy, 0, 2 * size, shift);
            Multiplex(proxy, 0, 0, 0, 2 * size, shift);
        }

        MostSignificantBit(proxy, 0, L_BIT * size, shift, false);

        for (int16_t i = L_BIT - 1; i >= 0; i--) {
            Compare(proxy, 0, 0, size, shift);
            Multiply(proxy, 0, 0, 2 * size, shift);
        }

        if(first_call) {
            Multiply(proxy, 0, 0, size, shift);
            Divide(proxy, 0, 0, size, shift, false);
            Multiplex(proxy, 0, 0, 0, size, shift);
        }
        return NULL;
    }
    return NULL;
}

/** Compute the division of a / b - not the integer approximation of the result
 *
 * @param proxy : Party instance
 * @param a : dividend
 * @param b : divider
 * @param first_call : indicates whether the aucDivide call is for the integer part of the division result, i.e. first call.
 * If it is the first call, then there will be the second call of aucDivide for the fractional part of the division result
 * @return a / b
 */
uint64_t Divide(Party *const proxy, uint64_t a, uint64_t b, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == helper) {
        Divide(proxy, nullptr, nullptr, 1, shift);
        return 0;
    } else {
        uint64_t* result_array = Divide(proxy, &a, &b, 1, shift);
        uint64_t result = result_array[0];
        delete[] result_array;
        return result;
    }
}

// TODO test
/** Perform division operation, or more specifically normalization operation, of two given inputs. The operation is
 * taken from SecureNN, but it is implemented by using the building blocks of CECILIA. Note that there is an implicit
 * assumption for Normalize to work correctly: the elements of a must be less than the corresponding elements of b.
 *
 * @param proxy
 * @param a: the nominators
 * @param b: the denominators
 * @param size: the number of elements in a and b
 * @return div: uint64_t vector consisting of elementwise division of a/b
 */
uint64_t* Normalise(Party *const proxy, const uint64_t *const a, const uint64_t *const b, size_t size, int shift = FRACTIONAL_BITS) {
    if (proxy->GetPRole() == proxy1 || proxy->GetPRole() == proxy2) {
        uint64_t *u = new uint64_t[size]; // holds how much needs to be subtracted from the nominator
        uint64_t *div = new uint64_t[size]; // resulting division
        for(size_t i = 0; i < size; i++) {
            u[i] = 0;
            div[i] = 0;
        }

        // iterate every bit of the fractional part to determine whether they are 1 or 0
        for(int i = 1; i <= shift; i++) {
            // compute the possible remaining of the nominator after subtracting denominator and previously subtracted value
            uint64_t *z = new uint64_t[size];
            for(size_t j = 0; j < size; j++) {
                z[j] = ((a[j] - u[j]) << i) - b[j];
            }

            uint64_t *msb_z = MostSignificantBit(proxy, z, size, shift);
            delete [] z;

            uint64_t *concat_cont_and_subt = new uint64_t[size * 2];
            uint64_t *twice_msb_z = new uint64_t[size * 2];
            for(size_t j = 0; j < size; j++) {
                twice_msb_z[j] = (proxy->GetPRole() << shift) - msb_z[j];
                twice_msb_z[j + size] = twice_msb_z[j];
                concat_cont_and_subt[j] = proxy->GetPRole() << (shift - i); // the contribution to the division result
                concat_cont_and_subt[j + size] = Truncate(proxy, b[j], i); // what to subtract from the nominator
            }
            delete [] msb_z;

            // computes possibly what to subtract and what to add & determines if we need to perform those operations
            uint64_t *tmp = Multiply(proxy, twice_msb_z, concat_cont_and_subt, 2 * size, shift);
            delete [] concat_cont_and_subt;
            delete [] twice_msb_z;

            for(size_t j = 0; j < size; j++) {
                div[j] = div[j] + tmp[j];
                u[j] = u[j] + tmp[j + size];
            }
            delete [] tmp;
        }

        delete [] u;
        return div;
    }
    else if (proxy->GetPRole() == helper) {
        for(int i = 1; i <= shift; i++) {
            MostSignificantBit(proxy, 0, size, shift);
            Multiply(proxy, 0, 0, 2 * size, shift);
        }
    }
    return nullptr;

}

#endif //CORE_H


