#include <iostream>
#include <queue>
#include <random>
#include <cmath>
#include <iomanip>

using namespace std;

// ********** Event Types **********

enum EventType {
    ARRIVAL,
    DEPARTURE
};

// ********** Event struct - each event has a time and a type **********

struct Event {
    double time;           // Time when event occurs
    EventType type;
    double arrival_time;   // זמן ההגעה של הלקוח

    // Constructor
    Event(double t, EventType tp, double at = 0.0) : time(t), type(tp), arrival_time(at) {}
    // comparison operator
    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// ********** MMN1Queue Class -  class for M/M/1/N queue **********

class MMN1Queue {
private:
    // Simulation parameters
    double lambda_rate;    // Arrival rate (λ)
    double mu_rate;        // Service rate (μ)
    int N;                 // Maximum capacity
    double T;              // Simulation time

    // current System state
    double current_time;
    int num_in_system;     // How many customers currently in system (queue + service)
    bool server_busy;      // Is the server busy

    // Statistics
    int customers_served;           // A - how many customers received full service
    int customers_blocked;          // Part of B - customers rejected because queue was full
    int customers_in_system_at_end; // Part of B - customers remaining at end
    double total_time_in_system;    // סכום זמני ההמתנה של כל הלקוחות

    // a priorety queue to handle the events , implemented with a vector of events
    // and the comparison func is the built in function greater
    // ( uses the comparison operator of the event => greater by time )
    //this will ensure that the event with the earliest time will be on the top
    priority_queue<Event, vector<Event>, greater<Event>> event_queue;

    // Random number generator - למשתנים אקראיים
    mt19937 rng;

public:

    // Constructor - Initialize the simulation

    MMN1Queue(double lambda, double mu, int n, double t, unsigned int seed = 0)
        : lambda_rate(lambda), mu_rate(mu), N(n), T(t),
          current_time(0.0), num_in_system(0), server_busy(false),
          customers_served(0), customers_blocked(0),
          customers_in_system_at_end(0), total_time_in_system(0.0) {

        // Initialize random number generator
        if (seed == 0) {
            random_device rd;
            rng.seed(rd());
        } else {
            rng.seed(seed);
        }
    }

 //In Poisson process with rate λ, inter-arrival times are exponentially
 //distributed with parameter λ. Using C++11's exponential_distribution
    /* from lectures : What is the distribution of the time interval between two Poisson arrivals?
    P(time between arrivals <= t)= 1 - P(time between arrivals> t)= 1 - P0(t) = 1 - exp(-λt)
             This is an exponential distribution
            */

    double generate_interarrival_time() {
        exponential_distribution<double> exp_dist(lambda_rate);
        return exp_dist(rng);
    }

    /*
     * Generate random service time - exponentially distributed with rate μ
     */
    double generate_service_time() {
        exponential_distribution<double> exp_dist(mu_rate);
        return exp_dist(rng);
    }

    //Schedule a new event in the event queue

    void schedule_event(double event_time, EventType event_type, double arrival_time = 0.0) {
        event_queue.push(Event(event_time, event_type, arrival_time));
    }

    /*
     * Handle arrival event
     */
    void handle_arrival() {
        if (num_in_system < N) {
            // There's room in the system
            num_in_system++;

            if (!server_busy) {
                // Server is idle - start service immediately
                server_busy = true;
                double service_time = generate_service_time();
                double departure_time = current_time + service_time;

                // Schedule the service completion
                schedule_event(departure_time, DEPARTURE, current_time);
            }
            // Otherwise - customer enters queue and waits
        } else {
            // No room - customer is blocked
            customers_blocked++;
        }
    }

    /*
     * Handle departure event
     */
    void handle_departure(Event& e) {
        num_in_system--;
        customers_served++;

        // חישוב זמן שהלקוח שהסתיים בילה במערכת
        total_time_in_system += current_time - e.arrival_time;

        // Check if there's another customer waiting
        if (num_in_system > 0) {
            // There are more customers - start serving the next one
            // The customer must have been waiting in queue, so start service now
            double service_time = generate_service_time();
            double departure_time = current_time + service_time;
            schedule_event(departure_time, DEPARTURE, current_time);
            // Server stays busy
        } else {
            // No more customers - server becomes idle
            server_busy = false;
        }
    }

    // the actual simulation process
    void run(int& A, int& B, double& avg_wait_time) {
        // Initialize: schedule the first arrival
        double first_arrival_time = generate_interarrival_time();
        schedule_event(first_arrival_time, ARRIVAL);

        // Main event loop
        while (!event_queue.empty()) {
            // Extract next event (earliest time)
            Event current_event = event_queue.top();
            event_queue.pop();

            // Check if we've passed simulation time
            if (current_event.time > T) {
                // Simulation time ended
                // Whatever remains in events won't happen
                break;
            }

            // Update current clock
            current_time = current_event.time;

            // Handle event based on type
            if (current_event.type == ARRIVAL) {
                handle_arrival();

                // Schedule next arrival (if we haven't reached T yet)
                double next_arrival_time = current_time + generate_interarrival_time();
                if (next_arrival_time <= T) {
                    schedule_event(next_arrival_time, ARRIVAL);
                }

            } else if (current_event.type == DEPARTURE) {
                handle_departure(current_event);
            }
        }

        // End of simulation - count who's left
        customers_in_system_at_end = num_in_system;

        // Calculate results
        A = customers_served;
        B = customers_blocked + customers_in_system_at_end;

        // חשב את זמן ההמתנה הממוצע
        avg_wait_time = (customers_served > 0) ? total_time_in_system / customers_served : 0.0;
    }
};


//for section 3.5
void run_experiment() {
    const double lambda = 10.0;
    const double mu = 15.0;
    const int N = 1000;
    const double theoretical_avg_wait = 1.0 / (mu - lambda); // M/M/1 approximation
    const double theoretical_A_per_T = lambda;               // E[A] = λ·T
    cout << fixed << setprecision(4);
    cout << "T   | Avg RelErr_A (%) | Avg RelErr_W (%)\n";
    cout << "-----------------------------------------\n";
    for (int T = 10; T <= 100; T += 10) {
        double sum_err_A = 0.0;
        double sum_err_W = 0.0;
        for (int run = 1; run <= 20; ++run) {
            MMN1Queue q(lambda, mu, N, T, run); // use run as seed
            int A, B;
            double avg_wait;
            q.run(A, B, avg_wait);
            double err_A = fabs(theoretical_A_per_T * T - A) /(theoretical_A_per_T * T) * 100.0;
            double err_W = fabs(theoretical_avg_wait - avg_wait) /theoretical_avg_wait * 100.0;
            sum_err_A += err_A;
            sum_err_W += err_W;
        }
        // הדפסת ממוצעים בלבד
        cout << setw(3) << T << " | "
             << setw(15) << sum_err_A / 20.0 << " | "
             << setw(15) << sum_err_W / 20.0 << "\n";
    }
}
//

// ********** Main function **********

int main(int argc, char* argv[]) {
//for section 3.5
    if (argc == 1) {
        run_experiment();
        return 0;
    }
//

    // Validate input
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " lambda mu N T" << endl;
        cerr << "  lambda: arrival rate (messages per time unit)" << endl;
        cerr << "  mu: service rate (messages per time unit)" << endl;
        cerr << "  N: number of messages in system (including one in service)" << endl;
        cerr << "  T: simulation run time" << endl;
        return 1;
    }

    // Read parameters
    double lambda_rate = atof(argv[1]);
    double mu_rate = atof(argv[2]);
    int N = atoi(argv[3]);
    double T = atof(argv[4]);

    // Validate values
    if (lambda_rate <= 0 || mu_rate <= 0 || N <= 0 || T <= 0) {
        cerr << "Error: All parameters must be positive" << endl;
        return 1;
    }

    // Create and run simulation
    // seed = 0 means we'll use random_device for a random seed
    MMN1Queue queue(lambda_rate, mu_rate, N, T, 0);
    int A, B;
    double avg_wait_time;
    queue.run(A, B, avg_wait_time);

    // Output according to requirements: A B + זמן ההמתנה הממוצע
    cout << fixed << setprecision(4);

    cout << "A = " << A << " B = " << B << " Avg_Waiting_Time = " << avg_wait_time << endl;

    return 0;
}
