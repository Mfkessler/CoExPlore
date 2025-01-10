# gunicorn_config_ranomics.py

# Bind to all interfaces on port 14569
bind = "0.0.0.0:14569"

# Number of worker processes
workers = 4

# Worker class (sync, gevent, eventlet, etc.)
worker_class = "sync"

# Maximum number of requests a worker will process before restarting
max_requests = 1000

# Timeout for requests
timeout = 30

# Logging
accesslog = "-"  # Access log to stdout
errorlog = "-"   # Error log to stdout
loglevel = "info"

# Daemonize the Gunicorn process (run in the background)
daemon = False
proc_name = "ranomics"
