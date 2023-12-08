import numpy as np
from .dna_language_specification.language import nucleotides
from .constraints.hairpin import Hairpin
from .constraints.constraints import Constraints
from .hyperparameters.hyperparameters import Hyperparameters

class ValidateKeyPayload:
    def __init__(self, constraints):
        self.payloads = set()
        self.motif_size = constraints.motif_size

        # Hairpin
        self.hairpin = Hairpin(constraints)

        # Homopolymer
        self.max_hom = constraints.max_hom
        self.hom_valid = True

        # GC-content
        self.min_gc = constraints.min_gc
        self.max_gc = constraints.max_gc

        # Keys
        self.keys = []
        self.key_size = constraints.key_size
        self.start_key_hom = []
        self.end_key_hom = []
        self.max_start_key_hom = {'A':0, 'T':0, 'C':0, 'G':0}
        self.max_end_key_hom = {'A':0, 'T':0, 'C':0, 'G':0}
        self.whole_key_hom_indices = {'A':-1, 'T':-1, 'C':-1, 'G':-1}
        self.min_key_gc_count = -1
        self.max_key_gc_count = -1

        # Payload
        self.payload_size = constraints.payload_size
        self.payload_num = constraints.payload_num
        self.max_start_payload_hom = {'A':0, 'T':0, 'C':0, 'G':0}
        self.max_end_payload_hom = {'A':0, 'T':0, 'C':0, 'G':0}
        self.max_payload_gc_count = -1
        self.min_payload_gc_count = -1

        # General Key Palyoad Size Validation
        self.valid_key_size = True
        self.valid_payload_size = True

    ### Get Validation ###

    def validate_all_constraints(self):
        if not (self.valid_payload_size and self.valid_key_size):
            return False
        hom_valid = self.get_homopolymer_validation()
        hairpin_valid = self.get_hairpin_validation()
        gc_valid = self.get_gc_validation()
        return hom_valid and hairpin_valid and gc_valid

    def validate_with_constraints(self, with_constraints):
        """This function validates a given set of payloads, list of key and set of 
        constraint thresholds for the constraints `with_constraints`.
        
        Parameters
        ----------
        with_constraints: set of str
            Set of strings containing a selection of the following constraints: 
            'hairpin', 'hom, 'gcContent'. Those will be the constraints that the
            palyoads and keys will have to conform to.
        
        Returns
        ----------
        valid: bool
            if `valid` is True, the set of keys and payloads confrom to the given constraints,
            else, if it is False, they violate at least one constraint.
        """
        if not (self.valid_payload_size and self.valid_key_size):
            return False
        for constraint in with_constraints:
            if constraint == 'hom':
                valid_hom = self.get_homopolymer_validation()
                if not valid_hom:
                    return False
            if constraint == 'gcContent':
                valid_gc = self.get_gc_validation()
                if not valid_gc:
                    return False
            if constraint == 'hairpin':
                valid_hairpin = self.get_hairpin_validation()
                if not valid_hairpin:
                    return False
        return True

    def get_homopolymer_validation(self):
        return self.valid_payload_size and self.valid_key_size and self.validate_hom()

    def get_hairpin_validation(self):
        return self.valid_payload_size and self.valid_key_size and self.validate_hairpin()

    def get_gc_validation(self):
        return self.valid_payload_size and self.valid_key_size and self.validate_gc()

    ##### Add Keys and Payloads #####

    ### Add Keys ###

    def set_keys(self, keys):
        self.keys = keys if isinstance(keys, list) else list(keys)
        successfully_added_keys = self.generate_key_pre_stats()
        self.valid_key_size = successfully_added_keys
        return successfully_added_keys
    
    ### Add Payloads ###
    
    def add_payload(self, new_payload):
        if len(new_payload) != self.payload_size:
            self.valid_payload_size = False
            return False
        self.payloads.add(new_payload)
        return self.add_hom_and_gc_stats(new_payload)
    
    def set_payloads(self, payloads):
        self.payloads = set()
        for payload in payloads:
            if not self.add_payload(payload):
                self.valid_payload_size = False
                return False
        return True

    ##### Pre-Stats #####

    ### Key Pre-Stats ###

    def generate_key_pre_stats(self):
        for j in range(len(self.keys)):
            if len(self.keys[j]) != self.key_size:
                return False
            key = self.keys[j]
            cur_base = key[0]
            cur_hom = 1
            cur_gc_count = 1 if cur_base in ['G', 'C'] else 0
            is_start = True

            for i in range(1, self.key_size):
                b = key[i]

                # Update start hom
                if cur_base != b:
                    if cur_hom > self.max_hom:
                        self.hom_valid = False
                    if is_start:
                        is_start = False
                        self.max_start_key_hom[cur_base] = max(self.max_start_key_hom[cur_base], \
                                                               cur_hom)
                        self.start_key_hom.append(cur_hom)
                    cur_base = b
                    cur_hom = 0
                
                cur_gc_count += 1 if b in ['G', 'C'] else 0

                cur_hom += 1
            
            # Update hom
            if cur_hom > self.max_hom:
                self.hom_valid = False
            if is_start:
                self.whole_key_hom_indices[cur_base] = j
                self.start_key_hom.append(cur_hom)
                self.max_start_key_hom[cur_base] = max(self.max_start_key_hom[cur_base], cur_hom)
            self.max_end_key_hom[cur_base] = max(self.max_end_key_hom[cur_base], cur_hom)
            self.end_key_hom.append(cur_hom)

            # Update overall gc count
            self.min_key_gc_count = cur_gc_count if self.min_key_gc_count == -1 \
                                                else min(self.min_key_gc_count, cur_gc_count)
            self.max_key_gc_count = max(self.max_key_gc_count, cur_gc_count)
        return True

    ### Payload Pre-Stats ###

    def add_hom_and_gc_stats(self, payload):
        is_start = True
        cur_base = payload[0]
        hom_len = 1
        cur_gc_count = 1 if cur_base in ['G', 'C'] else 0
        for i in range(1, self.payload_size):
            b = payload[i]

            # Update homopolymer stats
            if cur_base != b:
                if hom_len > self.max_hom:
                    self.hom_valid = False
                if is_start:
                    is_start = False
                    self.max_start_payload_hom[cur_base] = \
                                        max(self.max_start_payload_hom[cur_base], hom_len)
                hom_len = 0
                cur_base = b
            hom_len += 1

            cur_gc_count += 1 if b in ['G', 'C'] else 0

        # Update homopolymer stats
        if hom_len > self.max_hom:
            self.hom_valid = False
        if is_start:
            self.max_start_payload_hom[cur_base] = max(self.max_start_payload_hom[cur_base], \
                                                       hom_len)
        self.max_end_payload_hom[cur_base] = max(self.max_end_payload_hom[cur_base], hom_len)

        # Update GC stats
        self.max_payload_gc_count = max(self.max_payload_gc_count, cur_gc_count)
        self.min_payload_gc_count = cur_gc_count if self.min_payload_gc_count == -1 \
                                                else min(self.min_payload_gc_count, cur_gc_count)
        return True

    ##### Constraint Validations #####

    ### Homopolymer Validation ###

    def validate_hom(self):
        if len(self.payloads) == 0 or len(self.keys) == 0:
            return False
        if not self.hom_valid:
            return False
        for nuc in nucleotides:
            if not self.is_valid_homopolymer_for_base(nuc):
                return False
        return True

    def is_valid_homopolymer_for_base(self, base):
        # start and end keys
        if self.max_start_payload_hom[base] == self.payload_size:
            if self.whole_key_hom_indices[base] != -1:
                # whole motif hom with only one key
                if len(self.keys) == 1:
                    return False
                whole_key_index = self.whole_key_hom_indices[base]
                prev_key_index = (whole_key_index - 1) % len(self.keys)
                next_key_index = (whole_key_index + 1) % len(self.keys)
                end_prev_key_hom = self.end_key_hom[prev_key_index] \
                                        if self.keys[prev_key_index][self.key_size - 1] == base \
                                        else 0
                start_next_key_hom = self.start_key_hom[next_key_index] \
                                        if self.keys[next_key_index][0] == base else 0
                hom_len = end_prev_key_hom + 2 * self.key_size + start_next_key_hom +\
                          self.payload_size * 3
                if hom_len > self.max_hom:
                    return False
            else:
                max_prev_next_key = 0
                for i in range(len(self.keys)):
                    prev_key = 0
                    next_key = 0
                    if self.keys[i][self.key_size - 1] == base:
                        prev_key = self.end_key_hom[i] 
                    if self.keys[i][0] == base:
                        next_key = self.start_key_hom[i] 
                    if self.keys[(i + 1) % len(self.keys)][0] == base:
                        next_key = max(next_key, self.start_key_hom[(i + 1) % len(self.keys)])
                    max_prev_next_key = max(max_prev_next_key, prev_key + next_key)

                hom_len = max_prev_next_key + self.payload_size
                if hom_len > self.max_hom:
                    return False
        # start keys
        elif self.max_start_payload_hom[base] != 0:
            if self.whole_key_hom_indices[base] != -1:
                if self.max_end_payload_hom[base] == self.payload_size:
                    whole_key_index = self.whole_key_hom_indices[base]
                    prev_key_index = (whole_key_index - 1) % len(self.keys)
                    end_prev_key_hom = self.end_key_hom[prev_key_index] \
                                        if self.keys[prev_key_index][self.key_size - 1] == base \
                                        else 0
                    hom_len = self.max_start_payload_hom[base] + 2 * self.key_size + \
                              2 * self.payload_size + end_prev_key_hom
                    if hom_len > self.max_hom:
                        return False
                else:
                    hom_len = self.max_start_payload_hom[base] + self.key_size + \
                                self.max_end_payload_hom[base]
                    if hom_len > self.max_hom:
                        return False
            else:
                hom_len = self.max_start_payload_hom[base] + self.max_end_key_hom[base]
                if hom_len > self.max_hom:
                    return False
        # end keys
        elif self.max_end_payload_hom[base] != 0:
            if self.whole_key_hom_indices[base] != -1:
                if self.max_start_payload_hom[base] == self.payload_size:
                    whole_key_index = self.whole_key_hom_indices[base]
                    next_key_index = (whole_key_index + 1) % len(self.keys)
                    start_next_key_hom = self.start_key_hom[next_key_index] \
                                        if self.keys[next_key_index][0] == base else 0
                    hom_len = self.max_end_payload_hom[base] + 2 * self.key_size +\
                              2 * self.payload_size + start_next_key_hom
                    if hom_len > self.max_hom:
                        return False
                else:
                    hom_len = self.max_end_payload_hom[base] + self.key_size + \
                               self.max_start_payload_hom[base]
                    if hom_len > self.max_hom:
                        return False
            else:
                hom_len = self.max_end_payload_hom[base] + self.max_start_key_hom[base]
                if hom_len > self.max_hom:
                    return False
        return True

    ### GC Content Validation ###

    def validate_gc(self):
        if len(self.payloads) == 0 or len(self.keys) == 0:
            return False
        min_gc_count = np.ceil(self.min_gc * self.motif_size / 100)
        max_gc_count = np.floor(self.max_gc * self.motif_size / 100)
        if self.max_key_gc_count * 2 + self.max_payload_gc_count > max_gc_count or \
            self.min_key_gc_count * 2 + self.min_payload_gc_count < min_gc_count:
            return False
        return True

    ### Hairpin Validation ###

    def validate_hairpin(self):
        return self.hairpin.validate_hairpin(keys=self.keys, payloads=self.payloads)


if __name__ == '__main__':
    payload_size = 8
    payload_num = 5
    max_hom = 1
    max_hairpin = 2
    loop_size_min = 1
    loop_size_max = 1
    min_gc = 25
    max_gc = 60
    key_size = 1
    key_num = 5
    
    constraints = Constraints(payload_size=payload_size, payload_num=payload_num, \
                              max_hom=max_hom, max_hairpin=max_hairpin, \
                              min_gc=min_gc, max_gc=max_gc, key_size=key_size, \
                              key_num=key_num, loop_size_min=loop_size_min, \
                              loop_size_max=loop_size_max)

    with_constraints = {'hom', 'gcContent', 'hairpin'}

    keys = ['A', 'T']
    payloads = {'CGCGATCG', 'CGCTACTC', 'CGCTACAG', 'CATCGCAG', 'CGCTATGC'}
    print('keys: ', keys)
    print('payloads: ', payloads)
    validate = ValidateKeyPayload(constraints) 
    validate.set_keys(keys) 
    validate.set_payloads(payloads) 
    if validate.validate_with_constraints(with_constraints):
        print('Keys and payloads are valid!')
    else:
        print('Keys and payloads are not valid.')
