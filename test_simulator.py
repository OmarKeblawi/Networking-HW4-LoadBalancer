import subprocess
import time
import math
import sys

def get_mm1n_theory(lam, mu, Q, T):
    """
    Computes theoretical values for an M/M/1/N queue.
    N = Q + 1 (System capacity: queue + the one in service).
    """
    N = Q + 1
    rho = lam / mu
    ts_theory = 1 / mu

    # Case: Blocking Probability (Pn) [cite: 60, 62]
    if rho == 1:
        p_block = 1 / (N + 1)
        expected_L = N / 2
    else:
        p_block = (math.pow(rho, N) * (1 - rho)) / (1 - math.pow(rho, N + 1))
        expected_L = (rho / (1 - rho)) - ((N + 1) * math.pow(rho, N + 1) / (1 - math.pow(rho, N + 1)))

    # Effective arrival rate (lambda_eff) - arrivals that aren't dropped
    lam_eff = lam * (1 - p_block)
    
    expected_A = lam_eff * T
    expected_B = lam * T * p_block
    
    # Little's Law: W = L / lam_eff. Tw = W - Ts
    expected_W = expected_L / lam_eff if lam_eff > 0 else 0
    expected_Tw = max(0, expected_W - ts_theory)
    
    # Tend: After T, remaining items in system (expected_L) drain at rate mu
    expected_Tend = T + (expected_L / mu)

    return expected_A, expected_B, expected_Tend, expected_Tw, ts_theory

def run_test(name, T, M, P, lam, Q, mu):
    # Construct command: ./simulator T M P1..PM lam Q1..QM mu1..muM [cite: 104]
    cmd = ["./simulator", str(T), str(M)] + [str(p) for p in P] + [str(lam)] + \
          [str(q) for q in Q] + [str(m) for m in mu]
    
    print(f"--- [TEST] {name} ---")
    print(f"Command: {' '.join(cmd)}")
    start = time.time()
    try:
        # 2-minute limit 
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        elapsed = time.time() - start
    except subprocess.TimeoutExpired:
        print(f"FAILURE: Execution exceeded 2 minutes \n")
        return False

    if proc.returncode != 0:
        print(f"FAILURE: Simulator crashed or returned error {proc.returncode}\n")
        return False

    output = proc.stdout.strip().split()
    if len(output) != 5:
        print(f"FAILURE: Expected 5 output values, got {len(output)}\n")
        return False

    actual = [float(x) for x in output]
    
    # Calculate Theoretical (Summing independent M/M/1/N queues)
    t_A, t_B, t_Tend = 0, 0, 0
    t_Tw_numerator = 0
    t_Ts_numerator = 0
    total_lam_eff = 0

    for i in range(M):
        if P[i] > 0:
            # a=expected_A, b=expected_B, tend=expected_Tend, tw=expected_Tw, ts=ts_theory
            a, b, tend, tw, ts = get_mm1n_theory(lam * P[i], mu[i], Q[i], T)
            
            t_A += a
            t_B += b
            t_Tend = max(t_Tend, tend)
            
            # Calculate effective arrival rate for this server (lam * (1 - p_block))
            # This is roughly expected_A / T
            lam_eff_i = a / T
            
            # Weight averages by the number of messages actually served
            t_Tw_numerator += tw * lam_eff_i
            t_Ts_numerator += ts * lam_eff_i
            total_lam_eff += lam_eff_i

    # Final weighted theoretical averages
    final_t_Tw = t_Tw_numerator / total_lam_eff if total_lam_eff > 0 else 0
    final_t_Ts = t_Ts_numerator / total_lam_eff if total_lam_eff > 0 else 0

    metrics = ["A", "B", "Tend", "Tw", "Ts"]
    theory = [t_A, t_B, t_Tend, final_t_Tw, final_t_Ts]
    passed = True

    print(f"Time: {elapsed:.2f}s")
    for i in range(5):
        if metrics[i] == "B":
            # For B, test A + B is  within 10% margin
            passed_metric = abs((actual[0] + actual[i]) - (theory[0] + theory[i])) <= max((theory[0] + theory[i]) * 0.1, 0.05)
        else:
            # 10% margin 
            margin = max(theory[i] * 0.1, 0.05) if theory[i] != 0 else 0.05
            passed_metric = abs(actual[i] - theory[i]) <= margin
        status = "OK" if passed_metric else "FAIL"
        print(f"  {metrics[i]}: Act={actual[i]:.4f}, Exp={theory[i]:.4f} [{status}]")
        if status == "FAIL": passed = False
    
    print(f"RESULT: {'PASS' if passed else 'FAIL'}\n")
    return passed

# --- TEST VECTORS ---
test_vectors = [
    # 1. Sanity: Single Server Stable
    ("Single Server Stable", 5000, 1, [1.0], 20, [1000], [40]),
    
    # 2. Sanity: Zero Queue Blocking
    ("Zero Queue (M/M/1/1)", 5000, 1, [1.0], 20, [0], [40]),
    
    # 3. Sanity: Low Service Rate / Congested
    ("Congested Multi-Server", 5000, 4, [0.25]*4, 20, [100]*4, [0.5]*4),
    
    # 4. General: Balanced Load (λ=200, T=15000)
    ("Balanced High Load", 10000, 2, [0.5, 0.5], 200, [50, 50], [120, 120]),
    
    # 5. General: Asymmetric Probabilities
    ("Asymmetric Probabilities", 2000, 2, [0.8, 0.2], 100, [20, 20], [150, 50]),

    # 6. Edge Case: Single server receiving all traffic in M-server setup
    # Tests if LB correctly routes based on P=[1, 0, 0, 0]
    ("Multi-Server Static Routing", 5000, 4, [1.0, 0.0, 0.0, 0.0], 20, [1000]*4, [40]*4),

    # 7. Edge Case: Extremely high load with zero queue (Heavy Blocking) 
    # ρ = 2.5, N=1. High drop rate expected.
    ("Heavy Blocking (M/M/1/1)", 2000, 1, [1.0], 100, [0], [40]),

    # 8. Stability Threshold: λ is very close to μ
    # Tests system behavior near the stability limit (λ=14.9, μ=15.0).
    ("Stability Threshold", 500000, 1, [1.0], 14, [10000], [15.0]),

    # 9. Random Valid Case: Asymmetric Servers and Queues
    # Mixed capacities and service rates (λ=50)
    ("Asymmetric Mix", 3000, 3, [0.2, 0.5, 0.3], 50, [10, 50, 5], [20, 40, 15]),

    # 10. Edge Case: Probability distribution with very small values
    # Tests precision of routing for rare events.
    ("Micro Probabilities", 5000, 4, [0.001, 0.001, 0.997, 0.001], 20, [1000]*4, [40]*4),
    
    # 11. Edge Case: High Throughput/Short Duration
    # Testing performance with high λ (200) over a shorter T.
    ("High Velocity Burst", 500, 2, [0.6, 0.4], 200, [100, 100], [150, 100]),

    # 12. Edge Case: Extreme Asymmetry (One server very fast, one very slow)
    # Tests if Tw correctly averages across vastly different service times.
    ("Extreme Asymmetry", 50000, 2, [0.5, 0.5], 20, [100, 100], [100, 1]),

    # 13. Edge Case: Underutilized System
    # λ is very small compared to μ (ρ = 0.01). Tw should be near 0.
    ("Nearly Empty System", 5000, 1, [1.0], 1, [100], [100]),

    # 14. Edge Case: Very Short Time T
    # Tests if the simulator handles short durations without crashing.
    ("Short Duration", 10, 1, [1.0], 20, [100], [40]),

    # 15. Random High Capacity: Multiple servers, high λ, high μ
    ("High Capacity Cluster", 5000, 5, [0.2]*5, 500, [200]*5, [150]*5),

    # 16. Edge Case: High Queue, Low Lambda (Draining System)
    # Testing if Tend correctly handles a system that finishes almost exactly at T.
    ("Low Load High Buffer", 1000, 1, [1.0], 1, [5000], [100]),

    # 17. Random Valid: Many Servers, Mixed Probabilities
    # Tests the LB's ability to distribute load across 10 servers.
    ("Ten Server Array", 5000, 10, [0.1]*10, 200, [50]*10, [30]*10),

    # 18. Edge Case: Service Rate exactly equal to Arrival Rate (ρ=1)
    # Note: Expect high variance in Tw; this is the "unstable" boundary.
    ("Critical Stability (rho=1)", 500000, 1, [1.0], 20, [1000], [20]),

    # 19. Edge Case: Probability Zero for one server in a multi-server setup
    # Tests if the LB handles P=0 without division by zero or routing errors.
    ("Zero Probability Server", 5000, 2, [1.0, 0.0], 20, [100, 100], [40, 40]),

    # 20. Edge Case: Extremely Small Buffers (Q=1)
    # Different from Q=0; this tests if N=2 is handled correctly (1 in service, 1 in queue).
    ("Minimum Queue (Q=1)", 5000, 1, [1.0], 20, [1], [20]),

    # 21. Edge Case: Large T with Low Lambda (Efficiency Test)
    # Tests if the simulator can handle long idle periods without infinite loops.
    ("Sparse Arrivals", 10000, 1, [1.0], 0.05, [100], [10]),

    # 22. Edge Case: Saturated Capacity (All servers full)
    # All servers have mu=1, total lambda=10. High B and full queues at Tend.
    ("Total System Saturation", 2000, 3, [0.33, 0.33, 0.34], 10, [10, 10, 10], [1, 1, 1]),

    # 23. Edge Case: Non-Uniform mu and Q
    # Tests complex weighted averages for Tw and Ts.
    ("Complex Heterogeneous Servers", 5000, 2, [0.7, 0.3], 50, [20, 200], [60, 20])
]

if __name__ == "__main__":
    all_pass = True
    total_tests = len(test_vectors)
    passed_tests = 0
    failed_tests_names = []
    for vec in test_vectors:
        if run_test(*vec):
            passed_tests += 1
        else:
            all_pass = False
            failed_tests_names.append(vec[0])
    failed_tests = total_tests - passed_tests
    print(f"--- [SUMMARY] Tests Passed: {passed_tests}/{total_tests}, Failed: {failed_tests} ---")
    if failed_tests_names:
        print(f"Failed tests: {', '.join(failed_tests_names)}")
    sys.exit(0 if all_pass else 1)