#ifndef __GF_UTILS__
#define __GF_UTILS__

#include <vector>
using namespace std;

int sign(int n1) {
    if (n1 == 0) {
        return 1;
    }
    return n1 / abs(n1);
}

int maxElem(vector<int> &v) {
    int res = v[0];
    int size = v.size();
    for (int i = 0 ; i < size; i++) {
        res = max(res, v[i]);
    }
    return res;
}

int minElem(vector<int> &v) {
    int res = v[0];
    int size = v.size();
    for (int i = 0 ; i < size; i++) {
        res = min(res, v[i]);
    }
    return res;
}

#endif //__GF_UTILS__