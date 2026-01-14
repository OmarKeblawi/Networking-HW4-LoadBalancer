#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <iomanip>
#include <deque>

using namespace std;

// ********** Event Types **********
enum EventType {
    ARRIVAL,
    DEPARTURE
};

// ********** Event Struct **********
struct Event {
    double time;
    EventType type;
    int server_index;        // -1 for ARRIVAL
    double arrival_time;
    double service_duration;

    Event(double t, EventType tp, int s = -1, double a = 0.0, double sd = 0.0)
        : time(t), type(tp), server_index(s), arrival_time(a), service_duration(sd) {}

    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// ********** RNG **********
mt19937 rng;

// ********** Server **********
class Server {
public:
    int id;
    int capacity;           // Qi + 1
    double mu;

    int num_in_system;
    bool busy;
    deque<double> queue;

    long long served;
    long long dropped;
    double total_wait;
    double total_service;

    Server(int i, int Q, double mu_rate)
        : id(i), capacity(Q + 1), mu(mu_rate),
          num_in_system(0), busy(false),
          served(0), dropped(0),
          total_wait(0.0), total_service(0.0) {}

    double gen_service_time() {
        exponential_distribution<double> d(mu);
        return d(rng);
    }

    void accept(double t,
                priority_queue<Event, vector<Event>, greater<Event>>& pq) {
        if (num_in_system >= capacity) {
            dropped++;
            return;
        }

        num_in_system++;

        if (!busy) {
            busy = true;
            double s = gen_service_time();
            pq.push(Event(t + s, DEPARTURE, id, t, s));
        } else {
            queue.push_back(t);
        }
    }

    void depart(double t, double arr, double serv,
                priority_queue<Event, vector<Event>, greater<Event>>& pq,
                double T) {
        num_in_system--;

        // Do NOT count requests finishing after T
        if (t <= T) {
            served++;
            double start = t - serv;
            total_wait += (start - arr);
            total_service += serv;
        }

        if (!queue.empty()) {
            double next_arr = queue.front();
            queue.pop_front();
            double next_serv = gen_service_time();
            pq.push(Event(t + next_serv, DEPARTURE, id, next_arr, next_serv));
            busy = true;
        } else {
            busy = false;
        }
    }
};

// ********** Helpers **********
double gen_interarrival(double lambda) {
    exponential_distribution<double> d(lambda);
    return d(rng);
}

int pick_server(const vector<double>& P) {
    uniform_real_distribution<double> d(0.0, 1.0);
    double r = d(rng), acc = 0.0;
    for (size_t i = 0; i < P.size(); i++) {
        acc += P[i];
        if (r <= acc) return i;
    }
    return P.size() - 1;
}

// ********** Main **********
int main(int argc, char* argv[]) {
    random_device rd;
    rng.seed(rd());

    if (argc < 7) {
        cerr << "Usage: ./simulator T M P1..PM lambda Q1..QM mu1..muM\n";
        return 1;
    }

    int idx = 1;
    double T = atof(argv[idx++]);
    int M = atoi(argv[idx++]);

    if (argc != 4 + 3 * M) {
        cerr << "Error: wrong number of arguments for M=" << M << endl;
        return 1;
    }

    if (T <= 0 || M <= 0) {
        cerr << "Error: T > 0 and M > 0 required\n";
        return 1;
    }

    vector<double> P(M);
    double sumP = 0.0;
    for (int i = 0; i < M; i++) {
        P[i] = atof(argv[idx++]);
        if (P[i] < 0) {
            cerr << "Error: probabilities must be non-negative\n";
            return 1;
        }
        sumP += P[i];
    }

    if (fabs(sumP - 1.0) > 1e-6) {
        cerr << "Error: probabilities must sum to 1\n";
        return 1;
    }

    double lambda = atof(argv[idx++]);
    if (lambda < 0) {
        cerr << "Error: lambda must be >= 0\n";
        return 1;
    }

    vector<int> Q(M);
    for (int i = 0; i < M; i++) {
        Q[i] = atoi(argv[idx++]);
        if (Q[i] < 0) {
            cerr << "Error: Qi must be >= 0\n";
            return 1;
        }
    }

    vector<double> mu(M);
    for (int i = 0; i < M; i++) {
        mu[i] = atof(argv[idx++]);
        if (mu[i] <= 0) {
            cerr << "Error: mu_i must be > 0\n";
            return 1;
        }
    }

    vector<Server> servers;
    for (int i = 0; i < M; i++)
        servers.emplace_back(i, Q[i], mu[i]);

    priority_queue<Event, vector<Event>, greater<Event>> pq;

    double first = gen_interarrival(lambda);
    if (first <= T)
        pq.push(Event(first, ARRIVAL));

    double current = 0.0;
    double endT = 0.0;

    while (!pq.empty()) {
        Event e = pq.top();
        pq.pop();
        current = e.time;

        if (e.type == ARRIVAL) {
            if (current > T) continue;

            double next = current + gen_interarrival(lambda);
            if (next <= T)
                pq.push(Event(next, ARRIVAL));

            int s = pick_server(P);
            servers[s].accept(current, pq);
        } else {
            servers[e.server_index].depart(
                current, e.arrival_time, e.service_duration, pq, T);
            if (current > endT) endT = current;
        }
    }

    long long A = 0, B = 0;
    double Tw = 0.0, Ts = 0.0;

    for (auto& s : servers) {
        A += s.served;
        B += s.dropped;
        Tw += s.total_wait;
        Ts += s.total_service;
    }

    double avgTw = (A > 0) ? Tw / A : 0.0;
    double avgTs = (A > 0) ? Ts / A : 0.0;
    avgTw = max(0.0, avgTw);
    cout << fixed << setprecision(4);
    cout << A << " " << B << " " << endT << " "
         << avgTw << " " << avgTs << endl;

    return 0;
}
