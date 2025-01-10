import time
import os
import psycopg2

host = os.environ.get('DB_HOST', 'postgres')
user = os.environ.get('DB_USER', 'postgres')
password = os.environ.get('DB_PASSWORD', 'postgres')
dbname = os.environ.get('DB_NAME', 'wgcna')

while True:
    try:
        conn = psycopg2.connect(host=host, user=user, password=password, dbname=dbname)
        conn.close()
        print("Postgres is ready!")
        break
    except psycopg2.OperationalError:
        print("Waiting for Postgres...")
        time.sleep(5)
