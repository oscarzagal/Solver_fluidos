#include <bits/stdc++.h>

using namespace std;

int main() {

    vector<int> a = {23, 44, 44, 23, 11, 10, 1, 9, 11, 11};

    map<int, int> freq;

    for (const int x : a) {
        freq[x]++;

        cout << "freq[" << x << "] = " << freq[x] << "\n";
    }

    return 0;
}
