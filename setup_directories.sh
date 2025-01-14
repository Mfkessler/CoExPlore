#!/bin/bash

# Function to create directories with appropriate permissions
create_dir() {
  local dir="$1"
  echo "Creating directory: $dir"
  mkdir -p "$dir" || { echo "Error: Could not create $dir."; exit 1; }
  chmod 777 "$dir" || { echo "Error: Could not set permissions for $dir."; exit 1; }
}

# Create main directories
create_dir "output"
create_dir "app/app/static/images"
create_dir "app/output"

# Confirm successful creation
echo "All directories have been successfully created and permissions have been set."
