# BIF Template

Welcome to my BIF template, please follow the steps below. (Please get out)

## Project Description

Template for developing Python applications.

## Project Information

- **Name**: template thingy 
- **Description**: Simple python project template.
- **Version**: 0.0.1
- **Author**: Bjorn Wiggers (Bjornwiggers@hotmail.com)

## Setup Instructions

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

## Contact
   Bjorn Wiggers - bjornwiggers@hotmail.com or bjorn.wiggers@wur.nl
