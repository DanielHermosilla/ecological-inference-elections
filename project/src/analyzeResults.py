import json
import argparse
import numpy as np
import os
import re
from collections import defaultdict


def getGroupAndCandidate(filename):
    """
    Extracts (G, C) from filename in the format 'GXCXseedX.json'.
    """
    match = re.match(r"G(\d+)C(\d+)seed\d+\.json", filename)
    if match:
        return match.group(1), match.group(2)  # Returns (G, C) as strings
    return None, None


def readFilesAndWrite(input, output):

    if not os.path.exists(output):
        os.makedirs(output)  # Create output directory if it doesn't exist

    groupedData = defaultdict(list)

    # Read all JSON files and group by (G, C)
    for file in os.listdir(input):
        if file.endswith(".json"):
            G, C = getGroupAndCandidate(file)
            if G is None or C is None:
                print("An invalid file in the results directory was found")
                continue  # Skip invalid files
            inputPath = os.path.join(input, file)
            with open(inputPath, "r", encoding="utf-8") as f:
                data = json.load(f)
                groupedData[(G, C)].append(data)

                # Process each (G, C) group
    for (G, C), json_list in groupedData.items():
        mean_results = computeMeans(json_list)
        outputFilename = f"G{G}C{C}.json"
        outputPath = os.path.join(output, outputFilename)

        with open(outputPath, "w", encoding="utf-8") as f:
            json.dump(mean_results, f, indent=4)

        print(f"Processed Group {G}, Candidates {C} -> Saved to {outputPath}")


def computeMeans(json_list):
    """
    Computes mean for:
    - The last element of 'log_likelihood' across multiple JSONs in a group.
    - The 'time_taken' field.
    - The 'iterations_made' field.
    """
    logLikelihoodLastVal = []
    timeTakenLastVal = []
    iterationsMadeLastVal = []

    for data in json_list:
        try:
            # Extract last value of "log_likelihood"
            if "log_likelihood" in data and isinstance(data["log_likelihood"], list):
                last_value = data["log_likelihood"][-1]  # Extract last element
                if isinstance(last_value, (int, float)):
                    logLikelihoodLastVal.append(last_value)

            # Extract "time_taken" if it exists
            if "time_taken" in data and isinstance(data["time_taken"], (int, float)):
                timeTakenLastVal.append(data["time_taken"])

            # Extract "iterations_made" if it exists
            if "iterations_made" in data and isinstance(
                data["iterations_made"], (int, float)
            ):
                iterationsMadeLastVal.append(data["iterations_made"])
        except:
            print(f"There's an error in {data}")

    # Compute means
    mean_results = {
        "log_likelihood_mean_(last_iteration_value)": (
            np.mean(logLikelihoodLastVal) if logLikelihoodLastVal else None
        ),
        "time_taken_mean": np.mean(timeTakenLastVal) if timeTakenLastVal else None,
        "iterations_made_mean": (
            np.mean(iterationsMadeLastVal) if iterationsMadeLastVal else None
        ),
    }

    return mean_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute means from JSON files grouped by (G, C)."
    )
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input directory containing all of the results files (the one handed to ./util-exec)",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to the output directory where results will be saved (if it doesn't exist, it will be created)",
    )
    parser.add_argument(
        "method",
        type=str,
        help="The current methods are: `Multinomial`, `MVN PDF`, `MVN CDF`, `Hit and Run` and `Exact`",
    )

    if len(os.sys.argv) == 1:
        parser.print_help()
        print(
            'Remember that, when writing escaped arguments (arguments with spaces, such as `MVN PDF`) they should be surrounded by backticks "'
        )
        os.sys.exit(1)

    args = parser.parse_args()
    inputAppended = args.input + "/" + args.method + "/100B"
    outputAppended = args.output + "/" + args.method

    readFilesAndWrite(inputAppended, args.output)
