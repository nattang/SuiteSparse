import subprocess
import re
import json

def run_shell_command(command):
    process = subprocess.Popen(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    output_lines = []
    for line in process.stdout:
        print(line, end="")  # print live
        output_lines.append(line)

    process.wait()
    if process.returncode == 0:
        return "".join(output_lines).strip()
    else:
        print("Error: Command failed with return code", process.returncode)
        return None

# saves json of relevant stats
def parse_AMGX_output(output):
    # regexes for data
    obj_pattern = re.compile(r"Running CHOLMOD for ([\w\-]+)\.\.\.") 
    rows_pattern = re.compile(r"CHOLMOD sparse:\s+A:\s+(\d+)-by-\1, nz \d+,")
    analyze_time_pattern = re.compile(r"analyze\s+time:\s+(\d+\.\d+)\s*sec")
    solve_time_pattern = re.compile(r"solve\s+time:\s+([\d\.]+(?:e[+-]?\d+)?)\s*sec")
    error_pattern = re.compile(r"Total error between RXMesh and CHOLMOD:\s*([\d\.]+(?:e[+-]?\d+)?)")
    results = {}

    obj_names = obj_pattern.findall(output)
    rows = rows_pattern.findall(output)
    analyze_times = analyze_time_pattern.findall(output)
    solve_times = solve_time_pattern.findall(output)
    errors = error_pattern.findall(output)

    # print(obj_names, len(obj_names))
    # print(len(setup_times))

    for i, obj in enumerate(obj_names):
        num_rows = rows[i]
        analyze_time = float(analyze_times[i])
        solve_time = float(solve_times[i])
        error = float(errors[i])
        
        results[obj] = {
            "num_rows": num_rows,
            "total_analyze_time": analyze_time,
            "total_solve_time": solve_time,
            "total_error": error
        }

    return results



def main():
    command = "./cholmod_benchmark.sh" # run benchmarking solver
    output = run_shell_command(command)
    
    if output:
        print(output)
        results = parse_AMGX_output(output)
        output_path = "./results.json"
        with open(output_path, 'w') as json_file:
            json.dump(results, json_file, indent=4)
        print(json.dumps(results, indent=2))
    else:
        print("Failed to execute command.")

if __name__ == "__main__":
    main()
