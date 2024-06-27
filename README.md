# MEng Project

## About the Project

Currently, motif-based DNA storage is a cheap DNA storage option where only short DNA sequences are constructed and then assembled. However, due to biological and technological constraints involved in storing data in molecules, many encodings of arbitrary data to DNA end up being extremely error-prone, making DNA storage not a reliable storage option. Having to find an encoding design which takes into account those constraints as well as encodes a given data can be time extensive and challenging.

So to facilitate the search for a motif-based encoding design which conforms to a set of biological and technological constraints, we implemented:

1. A validation tool which verifies whether a set of motifs conforms to the given constraints.

The tool can be found on: http://ssb22.pythonanywhere.com

## Getting Started

### Prerequisites
* Python 3.8 or above.

### Installations
```bash
pip install -r /path/to/requirements.txt
```
## Usage
### Validation Tool

To verify whether a set of payloads and list of keys conform to a set of constraints, use the Validation Tool. 

The keys, payloads and constraints can be inputted directly in the main function of the validation_tool/key_payload_validation.py file.

To run the Validation Tool, run the following command from the root directory:
```bash
python -m validation_tool.key_payload_validation run
```

### Tests

To run the tests, first go to the unit_tests directory (validation_tool/unit_tests).

Run tests using the following command line: 
```bash
pytest filename.py
```
For example: 
```bash
pytest hairpin_tests.py
```
