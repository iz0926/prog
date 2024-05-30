#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<ll> vl;

mt19937 rng(random_device{}());
uniform_real_distribution<double> dist(0.0, 1.0);

ll A00(vl arr) {
    sort(arr.begin(), arr.end());
    while (arr.size() > 1) {
        ll tmp = arr[arr.size()-1] - arr[arr.size()-2];
        arr.pop_back(); arr.pop_back();
        arr.insert(lower_bound(arr.begin(), arr.end(), tmp), tmp);
    }
    return arr[0];
}

void A01(vl arr) {
    ll res = numeric_limits<ll>::max();
    for (int i = 0; i < 25000; i++) {
        ll tmp = 0;
        for (int j = 0; j < arr.size(); j++)
            dist(rng) < 0.5 ? tmp += arr[j] : tmp -= arr[j];
        res = min(res, abs(tmp));
    }
    cout << res << endl;
    return;
}

void A11(vl arr) {
    ll res = numeric_limits<ll>::max();
    for (int i = 0; i < 25000; i++) {
        vl val(arr.size(), 0);
        for (int j = 0; j < arr.size(); j++)
            val[int(dist(rng) * arr.size())] += arr[j];
        res = min(res, A00(val));
    }
    cout << res << endl;
    return;
}

void A02_03(vl arr, bool b) {
    int sta[arr.size()];
    ll curres = 0;
    for (int i = 0; i < arr.size(); i++) {
        sta[i] = dist(rng) < 0.5;
        curres += (2*sta[i]-1) * arr[i];
    }
    ll bestres = curres;
    for (int i = 0; i < 25000; i++) {
        int ind1 = int(dist(rng) * arr.size());
        int ind2 = int(dist(rng) * arr.size());
        if (ind1 == ind2) {
            ll tmp = curres + (2-4*sta[ind1]) * arr[ind1];
            sta[ind1] = 1-sta[ind1];
            double p = exp((abs(curres)-abs(tmp))/(1000000000*pow(0.8, i/300)));
            abs(tmp) < abs(curres) || (dist(rng) < p && b) ? curres = tmp : sta[ind1] = 1-sta[ind1]; 
        } else {
            ll tmp = curres + (2-4*sta[ind1]) * arr[ind1] + (2-4*sta[ind2]) * arr[ind2];
            sta[ind1] = 1-sta[ind1];
            sta[ind2] = 1-sta[ind2]; 
            double p = exp((abs(curres)-abs(tmp))/(1000000000*pow(0.8, i/300)));
            abs(tmp) < abs(curres) || (dist(rng) < p && b) ? curres = tmp : sta[ind1] = 1-sta[ind1], sta[ind2] = 1-sta[ind2]; 
        }
        abs(curres) < abs(bestres) ? bestres = curres : 0;
    }
    cout << abs(bestres) << endl;
    return;
}

void A12_13(vl arr, bool b) {
    int sta[arr.size()];
    vl val(arr.size(), 0);
    for (int i = 0; i < arr.size(); i++) {
        sta[i] = int(dist(rng) * arr.size());
        val[sta[i]] += arr[i];
    }
    ll res = A00(val), bestres = res;
    for (int i = 0; i < 25000; i++) {
        int ind1 = 0, ind2 = sta[ind1];
        while (sta[ind1] == ind2) {
            ind1 = int(dist(rng) * arr.size());
            ind2 = int(dist(rng) * arr.size());
        }
        // move ind1 to ind2
        int prev = sta[ind1];
        val[prev] -= 2*arr[ind1];
        val[ind2] += 2*arr[ind1];
        sta[ind1] = ind2;

        ll tmp = A00(val);
        double p = exp((res-tmp)/(1000000000*pow(0.8, i/300)));
        if (tmp <= res || (dist(rng) < p && b)) res = tmp;
        else { // Revert
            sta[ind1] = prev;
            val[prev] += 2*arr[ind1];
            val[ind2] -= 2*arr[ind1];
        }
        bestres = min(res, bestres);
    }
    cout << abs(bestres) << endl;
    return;
}

int main(int argc, char* argv[]) {
    if (argc != 4) return 0;
    string filename = argv[3];
    int alg = stoi(argv[2]);

    ifstream file(filename);
    string line;
    vl arr(100, 0);
    for (int i = 0; i < 100; i++) {
        getline(file, line);
        arr[i] = stoll(line);
    }

    if (alg == 0) cout << A00(arr) << endl;
    if (alg == 1) A01(arr);
    if (alg == 2) A02_03(arr, false);
    if (alg == 3) A02_03(arr, true);
    if (alg == 11) A11(arr);
    if (alg == 12) A12_13(arr, false);
    if (alg == 13) A12_13(arr, true);

    return 0;
}