from validation_tool.key_payload_validation import ValidateKeyPayload
from validation_tool.constraints.constraints import Constraints
import re

def validate_keys_and_payloads(keys_data, payloads_data, constraints, with_constraints):
    keys_data = keys_data.strip()
    payloads_data = payloads_data.strip()
    keys, key_size = check_keys(keys_data)
    if not keys:
        return "Wrong Input. Verify that all keys have the same size.", check_data(keys_data), check_data(payloads_data), False
    payloads, payload_size = check_payloads(payloads_data)
    if not payloads:
        return "Wrong Input. Verify that all payloads have the same size.", check_data(keys_data), check_data(payloads_data), False
    
    key_num = len(keys)
    payload_num = len(payloads)
    max_hom = constraints['maxHomopolymer']
    max_hairpin = constraints['maxHairpin']
    loop_size_min = constraints['loopMin']
    loop_size_max = constraints['loopMax']
    min_gc = constraints['gcContentMinPercentage']
    max_gc = constraints['gcContentMaxPercentage']

    if max_hairpin <= 0 or max_hom <= 0 or loop_size_min > loop_size_max or loop_size_min < 0 or \
        payload_size <= 0 or key_size <= 0 or key_num <= 0 or payload_num <= 0 or min_gc > max_gc \
        or min_gc < 0 or max_gc > 100:
        return "Wrong Input. Verify that constraints are correct.", check_data(keys_data), check_data(payloads_data), False

    constraints = Constraints(payload_size=payload_size, payload_num=payload_num, \
                             max_hom=max_hom, max_hairpin=max_hairpin, \
                             loop_size_min=loop_size_min, loop_size_max=loop_size_max,\
                             min_gc=min_gc, max_gc=max_gc, key_size=key_size, key_num=key_num)

    validate = ValidateKeyPayload(constraints) 
    validate.set_keys(keys) 
    validate.set_payloads(payloads) 
    isValid = validate.validate_with_constraints(with_constraints)
    for constraint in with_constraints:
        if constraint == "gcContent":
            isValid = validate.get_gc_validation()
            if not isValid:
                return "GC-Content constraints are violated.", check_data(keys_data), check_data(payloads_data), isValid
        elif constraint == "hom":
            isValid = validate.get_homopolymer_validation()
            if not isValid:
                return "Homopolymer constraints are violated.", check_data(keys_data), check_data(payloads_data), isValid
        elif constraint == "hairpin":
            isValid = validate.get_hairpin_validation()
            if not isValid:
                return "Hairpin constraints are violated.", check_data(keys_data), check_data(payloads_data), isValid
    return "", check_data(keys_data), check_data(payloads_data), True

def check_data(data):
    data = data.upper()
    data = re.sub(r"[^ATCG\s]+", "", data)
    data = re.sub(' +', ' ', data)
    return data

def check_payloads(data):
    """Removes all illegal characters of sequence."""

    data = data.upper()
    data = re.sub(r"[^ATCG ]+", "", data)
    data = re.sub(' +', ' ', data)
    data = data.strip()
    data = data.split(' ')
    payloads = set()
    payload_size = -1
    for payload in data:
        if len(payload) == 0 or (payload_size != -1 and payload_size != len(payload)):
            return False, False
        if payload_size == -1:
            payload_size = len(payload)
        payloads.add(payload)

    return payloads, payload_size

def check_keys(data):
    """Removes all illegal characters of sequence."""

    data = data.upper()
    data = re.sub(r"[^ATCG ]+", "", data)
    data = re.sub(' +', ' ', data)
    data = data.strip()
    data = data.split(' ')
    keys = []
    key_size = -1
    for key in data:
        if len(key) == 0 or (key_size != -1 and key_size != len(key)):
            return False, False
        if key_size == -1:
            key_size = len(key)
        keys.append(key)

    joints = []
    converse = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for k in keys:
        cur = k[0]
        joint = converse[k[0]]

        for i in range(1, len(k)):
            nuc = k[i]
            joint = converse[nuc] + joint
        joints.append(joint)
    return joints, key_size

