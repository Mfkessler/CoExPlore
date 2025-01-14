import os
import shutil
from sqlalchemy import create_engine, text
from sqlalchemy.exc import IntegrityError
import scanpy as sc
import ranomics.ortho as rortho
import ranomics.utils as rutils
from PIL import Image

def copy_figures_to_static(root_dir, target_dir='static/images'):
    """
    Copies all images from the 'figures' folders of the plant directories
    located in the specified root directory to the 'static/images/{PLANT}' target directory.
    
    Args:
    root_dir (str): Path to the root directory containing the plant directories.
    target_dir (str): Base path to the target directory for the copied images.
    """
    # Check if the root directory exists
    if not os.path.exists(root_dir):
        print(f"Error: The directory {root_dir} does not exist.")
        return
    
    # Create the target directory if it doesn't exist
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Iterate through all items in the root directory
    for item in os.listdir(root_dir):
        plant_dir = os.path.join(root_dir, item)
        
        # Check if suffix is "_ortho",
        # if os.path.isdir(plant_dir) and item.endswith("_ortho"):

        # Check if it is a directory with exactly two letters
        if os.path.isdir(plant_dir) and len(item) == 2:
            figures_dir = os.path.join(plant_dir, 'figures')
            
            # Check if a 'figures' folder exists
            if os.path.exists(figures_dir):
                dest_dir = os.path.join(target_dir, item)
                
                # Create the destination directory for the current plant
                if not os.path.exists(dest_dir):
                    os.makedirs(dest_dir)
                
                # Copy all files from the 'figures' folder
                for filename in os.listdir(figures_dir):
                    src_file = os.path.join(figures_dir, filename)
                    dest_file = os.path.join(dest_dir, filename)
                    if os.path.isfile(src_file):
                        shutil.copy(src_file, dest_file)
                        print(f"Copied: {src_file} to {dest_file}")
            else:
                print(f"No 'figures' folder found in {plant_dir}")


def copy_h5ad_to_data(root_dir, target_dir='data/h5ad'):
    """
    Copies all h5ad files from the 'analysis' subdirectories of plant directories and their corresponding ortho directories
    located in the specified root directory to 'data/h5ad' and 'data/h5ad/ortho' respectively.
    
    Args:
    root_dir (str): Path to the root directory containing the plant directories.
    target_dir (str): Base path to the target directory for the copied h5ad files.
    """
    # Ensure the root directory exists
    if not os.path.exists(root_dir):
        print(f"Error: The directory {root_dir} does not exist.")
        return
    
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Iterate through all items in the root directory
    for item in os.listdir(root_dir):
        plant_dir = os.path.join(root_dir, item)
        
        # Check if it is a directory and matches the pattern for plant or ortho folders
        if os.path.isdir(plant_dir):
            if len(item) == 2:  # Plant folder
                analysis_dir = os.path.join(plant_dir, 'analysis')

                # Dont't use analysis_dir
                h5ad_file_path = os.path.join(plant_dir, f"{item}.h5ad")
                if os.path.isfile(h5ad_file_path):
                    dest_file = os.path.join(target_dir, f"{item}.h5ad")
                    shutil.copy(h5ad_file_path, dest_file)
                    print(f"Copied: {h5ad_file_path} to {dest_file}")


def copy_tom_to_data(root_dir, target_dir='data/tom'):
    """
    Copies all tom files from the 'analysis' subdirectories of plant directories and their corresponding ortho directories
    located in the specified root directory to 'data/tom'.
    
    Args:
    root_dir (str): Path to the root directory containing the plant directories.
    target_dir (str): Base path to the target directory for the copied tom files.
    """
    # Ensure the root directory exists
    if not os.path.exists(root_dir):
        print(f"Error: The directory {root_dir} does not exist.")
        return
    
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Iterate through all items in the root directory
    for item in os.listdir(root_dir):
        plant_dir = os.path.join(root_dir, item)
        
        # Check if it is a directory and matches the pattern for plant or ortho folders
        if os.path.isdir(plant_dir):
            if len(item) == 2:  # Plant folder
                tom_file_path = os.path.join(plant_dir, f"tom_matrix.h5")
                if os.path.isfile(tom_file_path):
                    dest_file = os.path.join(target_dir, f"tom_matrix_{item}.h5")
                    shutil.copy(tom_file_path, dest_file)
                    print(f"Copied: {tom_file_path} to {dest_file}")


def generate_wgcna_db():
    """
    Generate a PostgreSQL database with the WGCNA browser data.
    """

    adatas = []
    names = ['AC', 'CS', 'EC', 'EG', 'HC', 'PS', 'SP', 'TT']

    for name in names:
        adata = sc.read_h5ad(f"/vol/share/ranomics_app/data/h5ad/{name}.h5ad")
        adatas.append(adata)

    df = rortho.transcript_ortho_browser("", adatas)

    engine = create_engine('postgresql+psycopg2://postgres:mokka1991@localhost/wgcna_data')

    df.to_sql('wgcna_browser', engine, if_exists='replace', index=True)


def generate_wgcna_db_cr():
    """
    Generate a PostgreSQL database with the WGCNA browser data for 'CR'.
    """

    adatas = []
    names = ['CR']

    for name in names:
        adata = sc.read_h5ad(f"/vol/share/ranomics_app/data/h5ad/{name}.h5ad")
        adatas.append(adata)

    # Generate the DataFrame
    df = rortho.transcript_ortho_browser("", adatas)

    # Check if 'ortho_id' column exists; if not, add it with empty strings
    if 'ortho_id' not in df.columns:
        df['ortho_id'] = ""
    if 'ortho_count' not in df.columns:
        df['ortho_count'] = 0

    # Save the DataFrame to PostgreSQL
    engine = create_engine('postgresql+psycopg2://postgres:mokka1991@localhost/wgcna_data')
    df.to_sql('wgcna_browser_cr', engine, if_exists='replace', index=True)

    print("Database wgcna_browser_cr created successfully!")


def generate_info_db():
    """
    Generate a PostgreSQL database with the WGCNA info data.
    """

    adatas = []
    names = ['AC', 'CS', 'EC', 'EG', 'HC', 'PS', 'SP', 'TT']

    for name in names:
        adata = sc.read_h5ad(f"/vol/share/ranomics_app/data/h5ad/{name}.h5ad")
        adatas.append(adata)

    engine = create_engine('postgresql+psycopg2://postgres:mokka1991@localhost/wgcna_data')
    # Retrieve info_df and sample_df
    info_df, sample_df = rutils.get_info_for_list(adatas)

    # Convert column names to lowercase and snake_case
    info_df.columns = [col.lower().replace(' ', '_') for col in info_df.columns]

    # Ensure that the 'species' column exists
    if 'species' not in info_df.columns:
        raise ValueError("The 'species' column is missing in info_df")

    # Drop the table if it exists
    with engine.connect() as conn:
        conn.execute(text("DROP TABLE IF EXISTS wgcna_info"))
        print("Dropped 'wgcna_info' table successfully (if it existed).")

    # Insert data into the new 'wgcna_info' table
    print("Creating and inserting 'wgcna_info' DataFrame into the database.")
    try:
        info_df.to_sql('wgcna_info', engine, if_exists='replace', index=False, method='multi')
        print("Inserted 'wgcna_info' DataFrame successfully.")
    except IntegrityError as e:
        print(f"IntegrityError while inserting 'wgcna_info' DataFrame: {e}")
    except Exception as e:
        print(f"An error occurred while inserting 'wgcna_info' DataFrame: {e}")


def generate_info_db_cr():
    """
    Generate a PostgreSQL database with the WGCNA info data for 'CR'.
    """

    adatas = []
    names = ['CR']

    for name in names:
        adata = sc.read_h5ad(f"/vol/share/ranomics_app/data/h5ad/{name}.h5ad")
        adatas.append(adata)

    engine = create_engine('postgresql+psycopg2://postgres:mokka1991@localhost/wgcna_data')

    info_df, _ = rutils.get_info_for_list(adatas)
    info_df.columns = [col.lower().replace(' ', '_') for col in info_df.columns]

    # Drop the table if it exists
    with engine.connect() as conn:
        conn.execute(text("DROP TABLE IF EXISTS wgcna_info_cr"))
        print("Dropped 'wgcna_info_cr' table successfully (if it existed).")

    info_df.to_sql('wgcna_info_cr', engine, if_exists='replace', index=False, method='multi')


def generate_dev_tables():
    """
    Generate PostgreSQL tables 'wgcna_info_all' and 'wgcna_browser_all' for all 9 plants.
    """
    # Liste aller Pflanzen
    names = ['AC', 'CS', 'EC', 'EG', 'HC', 'PS', 'SP', 'TT', 'CR']
    adatas = [sc.read_h5ad(f"/vol/share/ranomics_app/data/h5ad/{name}.h5ad") for name in names]

    # Verbindung zur PostgreSQL-Datenbank
    engine = create_engine('postgresql+psycopg2://postgres:mokka1991@localhost/wgcna_data')

    # Tabelle 'wgcna_browser_all' erstellen
    print("Creating 'wgcna_browser_all' table...")
    df_browser = rortho.transcript_ortho_browser("", adatas)
    if 'ortho_id' not in df_browser.columns:
        df_browser['ortho_id'] = ""
    if 'ortho_count' not in df_browser.columns:
        df_browser['ortho_count'] = 0

    df_browser.to_sql('wgcna_browser_all', engine, if_exists='replace', index=True, method='multi')
    print("Inserted 'wgcna_browser_all' table successfully.")

    # Tabelle 'wgcna_info_all' erstellen
    print("Creating 'wgcna_info_all' table...")
    info_df, _ = rutils.get_info_for_list(adatas)
    info_df.columns = [col.lower().replace(' ', '_') for col in info_df.columns]

    if 'species' not in info_df.columns:
        raise ValueError("The 'species' column is missing in info_df")

    # Alte Tabelle löschen, falls sie existiert
    with engine.connect() as conn:
        conn.execute(text("DROP TABLE IF EXISTS wgcna_info_all"))
        print("Dropped 'wgcna_info_all' table successfully (if it existed).")

    try:
        info_df.to_sql('wgcna_info_all', engine, if_exists='replace', index=False, method='multi')
        print("Inserted 'wgcna_info_all' table successfully.")
    except IntegrityError as e:
        print(f"IntegrityError while inserting 'wgcna_info_all': {e}")
    except Exception as e:
        print(f"An error occurred while inserting 'wgcna_info_all': {e}")


def update_prod_app(target_folder: str = "app_prod", source_folder: str = "app_dev"):
    """
    Update files in the target_folder with files from the source_folder.

    Parameters:
    - target_folder (str): The folder to update (default: "app_prod").
    - source_folder (str): The folder to copy files from (default: "app_dev").
    """
    # List of files to be copied and replaced
    files_to_update = [
        "cache.py",
        "factory.py",
        "tasks.py",
        "app.py",
        "celery_app.py",
        "extensions.py",
        "generate_index.py",
        "helper.py",
        "routes.py",
        os.path.join("static", "scripts.js"),
        os.path.join("static", "styles.css"),
    ]
    
    # Files in the 'templates' folder that should NOT be copied
    excluded_templates = {"index.html", "browser.html", "general_info.html"}
    
    # Check if the source directory exists
    if not os.path.exists(source_folder):
        print(f"Source folder '{source_folder}' does not exist.")
        return
    if not os.path.exists(target_folder):
        print(f"Target folder '{target_folder}' does not exist.")
        return
    
    # Copy and replace files
    for file in files_to_update:
        source_path = os.path.join(source_folder, file)
        target_path = os.path.join(target_folder, file)
        
        # Check if the source file exists
        if not os.path.exists(source_path):
            print(f"Source file '{source_path}' does not exist. Skipping...")
            continue
        
        # Create the target directory if it does not exist
        os.makedirs(os.path.dirname(target_path), exist_ok=True)
        
        # Copy the file
        try:
            shutil.copy2(source_path, target_path)
            print(f"Copied '{source_path}' to '{target_path}'.")
        except Exception as e:
            print(f"Error copying '{source_path}' to '{target_path}': {e}")
    
    # Copy files in the 'templates' folder
    templates_source = os.path.join(source_folder, "templates")
    templates_target = os.path.join(target_folder, "templates")
    
    if os.path.exists(templates_source):
        for template_file in os.listdir(templates_source):
            if template_file not in excluded_templates:
                source_path = os.path.join(templates_source, template_file)
                target_path = os.path.join(templates_target, template_file)
                
                # Copy the file
                try:
                    shutil.copy2(source_path, target_path)
                    print(f"Copied template '{source_path}' to '{target_path}'.")
                except Exception as e:
                    print(f"Error copying template '{source_path}' to '{target_path}': {e}")
    
    print(f"Update completed for target folder: '{target_folder}'.")


def generate_thumbnails(base_dir: str):
    """ 
    Generate thumbnails for all images in the specified directory.

    Parameters:
    - base_dir (str): The base directory containing the images.
    """
    
    if not os.path.exists(base_dir):
        print(f"Error: Directory '{base_dir}' does not exist.")
        return

    supported_formats = (".jpg", ".jpeg", ".png", ".bmp", ".tiff")

    for root, dirs, files in os.walk(base_dir):
        thumbnail_dir = os.path.join(root, 'thumbnails')
        os.makedirs(thumbnail_dir, exist_ok=True)

        for file in files:
            if file.lower().endswith(supported_formats):
                original_path = os.path.join(root, file)
                thumbnail_path = os.path.join(thumbnail_dir, file)

                if os.path.exists(thumbnail_path):
                    print(f"Thumbnail already exists: {thumbnail_path}")
                    continue

                try:
                    with Image.open(original_path) as img:
                        img.thumbnail((500, 500))
                        img.save(thumbnail_path, "PNG", quality=300)
                        print(f"Thumbnail created: {thumbnail_path}")
                except Exception as e:
                    print(f"Failed to create thumbnail for {original_path}: {e}")


def generate_all_thumbnails():
    """ 
    Generate thumbnails for all plant directories.
    """

    names = ['AC', 'CS', 'EC', 'EG', 'HC', 'PS', 'SP', 'TT', 'CR']
    dirs = [f"/vol/blast/ranomics-app/ranomics/flask/wgcna_app/app/static/images/{name}" for name in names]
    print("Generating thumbnails for the following directories:")
    print(dirs)

    for base_dir in dirs:
        print(f"Generating thumbnails in '{base_dir}'...")
        generate_thumbnails(base_dir)
        print("Thumbnail generation completed.")

if __name__ == "__main__":
    # copy_h5ad_to_data("../output/", target_dir="data/h5ad")
    # copy_figures_to_static("../output/", target_dir="static/images")
    # copy_tom_to_data("../output/", target_dir="data/tom")
    # generate_info_db()
    # generate_wgcna_db()
    # generate_info_db_cr()
    # generate_wgcna_db_cr()
    # generate_dev_tables()
    # update_prod_app(target_folder="app_cr", source_folder="app_dev")
    # update_prod_app(target_folder="app_prod", source_folder="app_dev")
    pass