# config.py


import json
import os


def is_dockerized() -> bool:
    """Check if the script is running inside a Docker container."""
    if os.path.exists("/.dockerenv"):
        return True

    if os.getenv("container") == "docker":
        return True

    try:
        with open("/proc/1/cgroup", "rt") as f:
            if "docker" in f.read() or "containerd" in f.read():
                return True
    except FileNotFoundError:
        pass

    try:
        with open("/proc/self/mounts", "rt") as f:
            if any("docker" in line or "containerd" in line for line in f):
                return True
    except FileNotFoundError:
        pass

    return False


def load_metadata_dict(metadata_json_path: str) -> dict:
    """
    Loads the metadata dictionary from a JSON file.

    Parameters:
    - metadata_json_path (str): Path to the metadata JSON file.

    Returns:
    - dict: The loaded metadata dictionary
    """

    try:
        with open(metadata_json_path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        print(
            f"Warning: {metadata_json_path} not found. Using an empty dictionary.")
        return {}
    except json.JSONDecodeError:
        print(
            f"Error: {metadata_json_path} contains invalid JSON. Using an empty dictionary.")
        return {}


# Load metadata based on environment
if is_dockerized():
    METADATA_DICT = load_metadata_dict("/metadata_dict.json")
else:
    METADATA_DICT = load_metadata_dict(
        "/vol/blast/wgcna/Project-Setup/wgcna-app/envs/metadata_dict.json")
