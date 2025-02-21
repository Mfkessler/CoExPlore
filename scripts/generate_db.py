import os
import json
import scanpy as sc
import wgcna.utils as rutils
import wgcna.ortho as rortho
from sqlalchemy import create_engine, text

def main():
    METADATA_DICT = json.loads(os.getenv("METADATA_DICT", "{}"))
    columns = list(METADATA_DICT.keys())
    h5ad_dir = os.environ.get("H5AD_DIR", "h5ad")

    # Environment variables for the DB from Docker Compose
    db_host = os.environ.get("DB_HOST", "localhost")
    db_user = os.environ.get("DB_USER", "postgres")
    db_password = os.environ.get("DB_PASSWORD", "postgres")
    db_name = os.environ.get("DB_NAME", "wgcna")

    # DB connection string
    engine = create_engine(f"postgresql+psycopg2://{db_user}:{db_password}@{db_host}/{db_name}")

    # Read AnnDatas in the h5ad directory
    h5ad_files = [f for f in os.listdir(h5ad_dir) if f.endswith(".h5ad")]
    adatas = [sc.read_h5ad(os.path.join(h5ad_dir, f)) for f in h5ad_files]

    # Create browser table
    df = rortho.transcript_ortho_browser("", adatas, possible_columns=columns)

    df.to_sql('wgcna_browser', engine, if_exists='replace', index=True)
    print("Table wgcna_browser created successfully!")

    # Create info table
    info_df, _ = rutils.get_info_for_list(adatas)
    info_df.columns = [col.lower().replace(' ', '_') for col in info_df.columns]

    with engine.connect() as conn:
        conn.execute(text("DROP TABLE IF EXISTS wgcna_info"))
        print("Dropped 'wgcna_info' table if it existed.")

    info_df.to_sql('wgcna_info', engine, if_exists='replace', index=False, method='multi')
    print("Table wgcna_info created successfully!")

    print("Database setup complete!")

if __name__ == "__main__":
    main()
