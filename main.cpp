#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <iomanip>
#include <deque>
#include <algorithm> // For max

using namespace std;

// ********** Event Types **********
enum EventType {
    ARRIVAL,    // Global arrival to Load Balancer
    DEPARTURE   // Departure from a specific server
};

// ********** Event Struct **********
struct Event {
    double time;
    EventType type;
    int server_index;       // -1 for ARRIVAL, 0 to M-1 for DEPARTURE
    double arrival_time;    // Timestamp when the request entered the system
    double service_duration;// Duration of the service (only relevant for DEPARTURE)

    // Constructor
    Event(double t, EventType tp, int s_idx = -1, double arr_t = 0.0, double serv_dur = 0.0)
        : time(t), type(tp), server_index(s_idx), arrival_time(arr_t), service_duration(serv_dur) {}

    // Priority Queue comparison (Min-Heap based on time)
    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// ********** Global Variables for Randomness **********
// We need a single source of randomness to ensure independent streams if needed,
// but for this simple simulation, a shared generator is sufficient.
mt19937 rng;

// ********** Server Class **********
class Server {
private:
    int id;
    int capacity;       // Q_i + 1 (Queue size + 1 for service)
    double mu_rate;     // Service rate
    
    // Internal state
    int num_in_system;  // Current number of requests (Queue + Service)
    bool busy;          // Is the server currently processing?
    deque<double> waiting_queue; // Stores arrival times of requests in queue

public:
    // Statistics for this server
    int served_count;
    int dropped_count; // Dropped due to full queue
    double total_wait_time; // Sum of (Start Service Time - Arrival Time)
    double total_service_time; // Sum of Service Durations

    Server(int server_id, int Q, double mu) 
        : id(server_id), capacity(Q + 1), mu_rate(mu), 
          num_in_system(0), busy(false),
          served_count(0), dropped_count(0), total_wait_time(0), total_service_time(0) {}

    // Attempt to accept a request
    // Returns true if accepted, false if dropped (queue full)
    bool accept_arrival(double current_time, priority_queue<Event, vector<Event>, greater<Event>>& pq) {
        if (num_in_system >= capacity) {
            dropped_count++;
            return false;
        }

        num_in_system++;

        if (!busy) {
            // Start service immediately
            busy = true;
            double service_time = generate_service_time();
            
            // Schedule departure
            // Wait time is 0 because it started immediately
            // We store service_time in the event to log it later
            pq.push(Event(current_time + service_time, DEPARTURE, id, current_time, service_time));
        } else {
            // Queue the request
            waiting_queue.push_back(current_time);
        }
        return true;
    }

    // Handle departure
    void handle_departure(double current_time, double arrival_time, double service_duration, 
                          priority_queue<Event, vector<Event>, greater<Event>>& pq) {
        num_in_system--;
        served_count++;
        
        // Update stats
        // Wait Time = (Departure Time - Service Duration) - Arrival Time
        // Effectively: Start_Service_Time - Arrival_Time
        double start_service_time = current_time - service_duration;
        double wait_time = start_service_time - arrival_time;
        
        total_wait_time += wait_time;
        total_service_time += service_duration;

        // Check if there are more requests in queue
        if (!waiting_queue.empty()) {
            // Start processing next one
            double next_arrival_time = waiting_queue.front();
            waiting_queue.pop_front();
            
            double next_service_time = generate_service_time();
            
            // Schedule next departure
            pq.push(Event(current_time + next_service_time, DEPARTURE, id, next_arrival_time, next_service_time));
            // busy remains true
        } else {
            busy = false;
        }
    }

    double generate_service_time() {
        exponential_distribution<double> exp_dist(mu_rate);
        return exp_dist(rng);
    }
};

// ********** Helper Functions **********

double generate_interarrival_time(double lambda) {
    exponential_distribution<double> exp_dist(lambda);
    return exp_dist(rng);
}

int pick_server(const vector<double>& probs) {
    uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(rng);
    double cumulative = 0.0;
    for (size_t i = 0; i < probs.size(); ++i) {
        cumulative += probs[i];
        if (r <= cumulative) {
            return i;
        }
    }
    return probs.size() - 1; // Fallback for rounding errors
}

// ********** Main **********

int main(int argc, char* argv[]) {
    // Seed RNG
    random_device rd;
    rng.seed(rd());

    // 1. Parse Arguments
    // Expected: ./simulator T M P1...PM lambda Q1...QM mu1...muM
    
    // Minimal check: T(1) + M(1) + P(1) + lambda(1) + Q(1) + mu(1) + prog(1) = 7
    if (argc < 7) {
        cerr << "Usage: ./simulator T M P1...PM lambda Q1...QM mu1...muM" << endl;
        return 1;
    }

    int arg_idx = 1;
    double T = atof(argv[arg_idx++]);
    int M = atoi(argv[arg_idx++]);

    // Check if we have enough arguments based on M
    // We need M probs, 1 lambda, M Qs, M mus
    if (argc != (3 + M + 1 + M + M)) {
        cerr << "Error: Incorrect number of arguments for M=" << M << endl;
        return 1;
    }

    vector<double> probs;
    double prob_sum = 0.0;
    for (int i = 0; i < M; ++i) {
        probs.push_back(atof(argv[arg_idx++]));
        prob_sum += probs.back();
    }
    if(abs(prob_sum - 1.0) > 1e-6) {
        cerr << "Error: Probabilities must sum to 1." << endl;
        return 1;
    }

    double lambda = atof(argv[arg_idx++]);

    vector<int> queues;
    for (int i = 0; i < M; ++i) queues.push_back(atoi(argv[arg_idx++]));

    vector<double> mus;
    for (int i = 0; i < M; ++i) mus.push_back(atof(argv[arg_idx++]));

    // 2. Setup System
    vector<Server> servers;
    for (int i = 0; i < M; ++i) {
        servers.emplace_back(i, queues[i], mus[i]);
    }

    priority_queue<Event, vector<Event>, greater<Event>> pq;

    // 3. Simulation Loop
    double current_time = 0.0;
    double tend = 0.0; // Time of last departure

    // Schedule first arrival
    double first_arr = generate_interarrival_time(lambda);
    if (first_arr <= T) {
        pq.push(Event(first_arr, ARRIVAL));
    }

    while (!pq.empty()) {
        Event e = pq.top();
        pq.pop();

        // Update time
        current_time = e.time;

        if (e.type == ARRIVAL) {
            // ARRIVALS stop after time T
            if (current_time > T) continue;

            // 1. Schedule next arrival
            double next_arr = current_time + generate_interarrival_time(lambda);
            if (next_arr <= T) {
                pq.push(Event(next_arr, ARRIVAL));
            }

            // 2. Route request to a server
            int server_idx = pick_server(probs);
            
            // 3. Try to add to server
            servers[server_idx].accept_arrival(current_time, pq);

        } else if (e.type == DEPARTURE) {
            // Handle departure at specific server
            int s_idx = e.server_index;
            servers[s_idx].handle_departure(current_time, e.arrival_time, e.service_duration, pq);
            
            // Track last departure time
            if (current_time > tend) tend = current_time;
        }
    }

    // 4. Aggregate Results
    long long total_A = 0;
    long long total_B = 0;
    double total_wait_sum = 0.0;
    double total_service_sum = 0.0;

    for (const auto& s : servers) {
        total_A += s.served_count;
        total_B += s.dropped_count;
        total_wait_sum += s.total_wait_time;
        total_service_sum += s.total_service_time;
    }

    double avg_Tw = (total_A > 0) ? (total_wait_sum / total_A) : 0.0;
    double avg_Ts = (total_A > 0) ? (total_service_sum / total_A) : 0.0;

    // 5. Output
    cout << fixed << setprecision(4);
    cout << total_A << " " << total_B << " " << tend << " " << avg_Tw << " " << avg_Ts << endl;

    return 0;
}