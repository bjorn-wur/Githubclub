# BIF Template

Welcome to BIF template, please follow the steps below.

## Project Description

Template for developing Python applications.

## Project Information

- **Name**: template thingy 
- **Description**: Simple python project template.
- **Version**: 0.0.1
- **Author**: Bjorn Wiggers

## Quickstart guide:

   ```bash
   git clone git@github.com:bjorn-wur/Githubclub.git

   python -m venv env
   source env/bin/activate  # On Windows, use `env\Scripts\activate`
   pip install -r requirements.txt
   pip install -r dev-requirements.txt

   pre-commit install
   ```

if you already have a repo but are unable to push please change the remote origin via the following command.
git remote set-url origin git@git.wur.nl:wigge031/Githubclub.git

## Full Setup Instructions:

Follow these steps to set up the project on your local machine:

1. **Clone the Repository**

   Start by cloning the project repository (Not Applicable currently):

   ```bash
   git clone <repository_url>
   ```

   or just copy it i dont care.

2. **Set Up a Virtual Environment**

   It is recommended to use a virtual environment to manage project dependencies. You can create and activate a virtual environment using:

   ```bash
   python -m venv env
   source env/bin/activate  # On Windows, use `env\Scripts\activate`
   ```

3. **Install Dependencies**

   Once the virtual environment is activated, you can install the required dependencies. There are two types of dependencies: main requirements and development requirements.

   - **Main Requirements**: Install the main project dependencies listed in `requirements.txt`:

     ```bash
     pip install -r requirements.txt
     ```

   - **Development Requirements**: Install the development dependencies listed in `dev-requirements.txt`:
      currently there are none

     ```bash
      pip install -r dev-requirements.txt
     ```

4. **Set Up Pre-commit Hooks**

   The project uses `pre-commit` to maintain code quality. Set up the pre-commit hooks using:

   ```bash
      pre-commit install
   ```

5. **Run Tests**

   To ensure everything is set up correctly, you can run unit tests in .\test directory

   ```bash
   python -m pytest
   ```


6. **requirements.txt**
   to keep thing simple just use requirements freeze in order to add requirements to the text file
   (so pip freeze > requirements.txt)


## Additional Information
you might have to add credentials:
   git config --global user.email bjorn.wiggers@wur.nl
   git config --global user.name bjorn-wur

pre-commit hooks can be run without commit with pre-commit run --all-files
all local modules in the enviroment can be saved using. 
pip freeze > requirements.txt

## Contact
   Bjorn Wiggers - bjorn.wiggers@wur.nl
