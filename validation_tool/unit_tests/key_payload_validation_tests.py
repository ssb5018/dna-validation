import pytest
import numpy as np
from .. import key_payload_validation as v
from ..constraints import constraints as c
from ..hyperparameters import hyperparameters as h

def get_constraints(maxHairpin=1, loopSize=-1, payloadSize=10, keySize=2, payloadNum=5, maxHom=1, minGC=25, maxGC=60, loopSizeMin=1, loopSizeMax=1):
    payloadSize = payloadSize
    payloadNum = payloadNum
    maxHom = maxHom
    maxHairpin = maxHairpin
    loopSize = loopSize
    loopSizeMin = loopSizeMin
    loopSizeMax = loopSizeMax
    minGc = minGC
    maxGc = maxGC
    keySize = keySize
    keyNum = 10 
    
    return c.Constraints(payloadSize, payloadNum, maxHom, maxHairpin, loopSize, minGc, maxGc, keySize, keyNum, loopSizeMin, loopSizeMax)

###### Homopolymer Validation ######

def test_valid_homopolymer_empty_returns_false():
    maxHom = 1
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = []
    payloads = set()
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = False
    assert result == valid

def test_valid_homopolymer_wrong_length_returns_false():
    maxHom = 1
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['ATC']
    payloads = {'TATGCTGCTG'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = False
    assert result == valid

def test_valid_homopolymer():
    maxHom = 1
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['AG']
    payloads = {'TATGCTGCTG'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = True
    assert result == valid

def test_valid_homopolymer_false_within_key():
    maxHom = 1
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['AA']
    payloads = {'TATGCTGCTG'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = False
    assert result == valid

def test_valid_homopolymer_false_within_payload():
    maxHom = 1
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['AG']
    payloads = {'TATCCTGCTG'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = False
    assert result == valid

def test_valid_homopolymer_false_between_key_payload():
    maxHom = 1
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['AT']
    payloads = {'TATGCTGCTG'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = False
    assert result == valid

def test_valid_homopolymer_true_multiple_keys_order_matters():
    maxHom = 11
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['AC', 'GA', 'CA']
    payloads = {'CCCCCCCCCC'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = True
    assert result == valid

def test_valid_homopolymer_true_multiple_keys_false():
    maxHom = 11
    constraints = get_constraints(maxHom=maxHom)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CA', 'GA', 'AC']
    payloads = {'CCCCCCCCCC'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_homopolymer_validation()
    result = False
    assert result == valid

###### GC-Content Validation ######

def test_valid_gc():
    minGC = 25
    maxGC = 58
    constraints = get_constraints(minGC=minGC, maxGC=maxGC)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CA', 'GA', 'AC']
    payloads = {'CCACACACAC'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_gc_validation()
    result = True
    assert result == valid

def test_valid_gc_false():
    minGC = 25
    maxGC = 58
    constraints = get_constraints(minGC=minGC, maxGC=maxGC)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CA', 'GG', 'AC']
    payloads = {'CCACACACAC'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_gc_validation()
    result = False
    assert result == valid

###### Hairpin Validation ######

def test_valid_hairpin():
    constraints = get_constraints()
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CA', 'GG', 'AC']
    payloads = {'CCACACACAC'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = True
    assert result == valid

def test_valid_hairpin_false_key_payload():
    constraints = get_constraints()
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CA', 'GG', 'AC']
    payloads = {'CGACACACAC'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = False
    assert result == valid

def test_valid_hairpin_false_in_payload():
    constraints = get_constraints()
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CA', 'GG', 'AC']
    payloads = {'CCACGCAGCA'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = False
    assert result == valid

def test_valid_hairpin_false_in_key():
    keySize = 5
    constraints = get_constraints(keySize=keySize)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CGACG', 'AAACA', 'AAAAA']
    payloads = {'CCACACAGCA'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = False
    assert result == valid

def test_valid_hairpin_false_same_key():
    keySize = 5
    loopSizeMin = 13
    loopSizeMax = 13
    constraints = get_constraints(keySize=keySize, loopSizeMin=loopSizeMin, loopSizeMax=loopSizeMax)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CGACG', 'AACAA', 'AAAAA']
    payloads = {'CCACACAGCA'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = False
    assert result == valid

def test_valid_hairpin_false_order_key_matters():
    keySize = 5
    loopSizeMin = 44
    loopSizeMax = 44
    constraints = get_constraints(keySize=keySize, loopSizeMin=loopSizeMin, loopSizeMax=loopSizeMax)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CGAAG', 'AACAA', 'ACGAA']
    payloads = {'CCACACAGCA'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = False
    assert result == valid

def test_valid_hairpin_false_in_different_payloads():
    keySize = 5
    loopSizeMin = 12
    loopSizeMax = 12
    constraints = get_constraints(keySize=keySize, loopSizeMin=loopSizeMin, loopSizeMax=loopSizeMax)
    validation = v.ValidateKeyPayload(constraints)
    keys = ['CGAAG', 'AACAA', 'ACGAA']
    payloads = {'CCACACAGCA', 'CCCGACAGCA', 'CCACGCAGCA'}
    validation.set_keys(keys)
    validation.set_payloads(payloads)
    valid = validation.get_hairpin_validation()
    result = False
    assert result == valid
