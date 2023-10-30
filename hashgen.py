import hashlib

# Global dictionary to store the mappings between unique IDs and strings
id_to_string_mapping = {}

def generate_unique_id(input_string, length=3):
    # Use a hash function (e.g., MD5) to generate a hash from the input string.
    # You can choose a different hash function if needed.
    hash_object = hashlib.md5(input_string.encode())
    unique_id = hash_object.hexdigest()[:length]
    
    # Store the mapping between the unique ID and the original string
    id_to_string_mapping[unique_id] = input_string
    
    return unique_id

def convert_id_to_string(unique_id):
    # Look up the original string based on the unique ID
    return id_to_string_mapping.get(unique_id, None)

# Test the functions
#input_string = "12345678babb"
#unique_id = generate_unique_id(input_string, length=3)
#
## Convert the ID back to the original string
#original_string = convert_id_to_string(unique_id)