#ifndef LICOM3_KOKKOS_SRC_CPP_IMPL_CPP_HALO_TRANSPOSE_HPP_
#define LICOM3_KOKKOS_SRC_CPP_IMPL_CPP_HALO_TRANSPOSE_HPP_

template<typename T>
void cpp_get_halo_transpose (const T* arrSrc, T* const arrObj,
    const int &startLayer, const int &lenLayer, 
        const int &lenA, const int &lenB, const int &lenC) {

  int startB, endB, startC, endC;

  for (int a = 0; a < lenA; ++a) {

    startB = startLayer;
    endB   = lenB - startLayer;
    startC = startLayer;
    endC   = startLayer + lenLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
      }
    }

    startB = startLayer;
    endB   = startLayer + lenLayer;
    startC = startLayer + lenLayer;
    endC   = lenC - startLayer - lenLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
      }
    }

    startB = lenB - startLayer - lenLayer;
    endB   = lenB - startLayer;
    startC = startLayer + lenLayer;
    endC   = lenC - startLayer - lenLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
      }
    }

    startB = startLayer;
    endB   = lenB - startLayer;
    startC = lenC - startLayer - lenLayer;
    endC   = lenC - startLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
      }
    }
  }

  return ;
}

template<typename T>
void cpp_put_halo_transpose (const T* arrSrc, T* const arrObj,
    const int &startLayer, const int &lenLayer, 
        const int &lenA, const int &lenB, const int &lenC) {

  int startB, endB, startC, endC;
  for (int a = 0; a < lenA; ++a) {

    startB = startLayer;
    endB   = lenB - startLayer;
    startC = startLayer;
    endC   = startLayer + lenLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[a*lenB*lenC + b*lenC + c] = arrSrc[b*lenC*lenA + c*lenA + a];
      }
    }

    startB = startLayer;
    endB   = startLayer + lenLayer;
    startC = startLayer + lenLayer;
    endC   = lenC - startLayer - lenLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[a*lenB*lenC + b*lenC + c] = arrSrc[b*lenC*lenA + c*lenA + a];
      }
    }

    startB = lenB - startLayer - lenLayer;
    endB   = lenB - startLayer;
    startC = startLayer + lenLayer;
    endC   = lenC - startLayer - lenLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[a*lenB*lenC + b*lenC + c] = arrSrc[b*lenC*lenA + c*lenA + a];
      }
    }

    startB = startLayer;
    endB   = lenB - startLayer;
    startC = lenC - startLayer - lenLayer;
    endC   = lenC - startLayer;
    for (int b = startB; b < endB; ++b) {
      for (int c = startC; c < endC; ++c) {
        arrObj[a*lenB*lenC + b*lenC + c] = arrSrc[b*lenC*lenA + c*lenA + a];
      }
    }
  }

  return ;
}
#endif // LICOM3_KOKKOS_SRC_CPP_IMPL_CPP_HALO_TRANSPOSE_HPP_