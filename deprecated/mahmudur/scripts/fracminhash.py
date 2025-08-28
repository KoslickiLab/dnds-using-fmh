class FracMinHash:
    '''
    FracMinHash class.
    '''
    def __init__(self, scale_factor, max_hash_value, initial_set=None):
        '''
        Create an FMH with given scale_factor in (0,1), largest hash value H.
        If initial_set is provided, that needs to be constructed as a set(),
        and must only contain integers in [0,H)
        Returns: None
        '''
        if initial_set is None:
            self.hash_set = set()
        else:
            self.hash_set = set(initial_set)
        self.H = max_hash_value
        self.scale_factor = scale_factor

    def add_value(self, hash_value):
        '''
        Add a hash value to the sketch.
        Returns: None
        '''
        if hash_value <= self.H * self.scale_factor:
            self.hash_set.add(hash_value)
        
    def add_values(self, hash_values):
        '''
        Add multiple hash values to the sketch.
        Returns: None
        '''
        for hash_value in hash_values:
            self.add_value(hash_value)

    def remove(self, hash_value):
        '''
        Remove a hash value from sketch.
        Returns: None
        '''
        self.hash_set -= hash_value

    def get_containment(self, smh):
        '''
        Obtain containment of provided sketch within this sketch.
        Returns: float -- containment value
        '''
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / len(self.hash_set)

    def get_scaled_containment(self, smh, num_all_elements):
        '''
        Obtain scaled containment of provided sketch within this sketch.
        Accounts for the bias factor. Must provide L from the paper.
        Returns: float -- containment value
        '''
        bf = 1.0 - (1.0 - self.scale_factor) ** int(num_all_elements)
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / ( len(self.hash_set) * bf )

    def get_sketch_size(self):
        '''
        Returns: int -- size of the sketch
        '''
        return len( self.hash_set )
