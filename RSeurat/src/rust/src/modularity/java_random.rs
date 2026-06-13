/// Java-compatible linear congruential generator (java.util.Random semantics).
#[derive(Clone)]
pub struct JavaRandom {
    seed: u64,
}

impl JavaRandom {
    pub fn new(seed: u64) -> Self {
        let mut rng = Self { seed: 0 };
        rng.set_seed(seed);
        rng
    }

    pub fn set_seed(&mut self, seed: u64) {
        self.seed = (seed ^ 0x5DEECE66D) & ((1u64 << 48) - 1);
    }

    fn next(&mut self, bits: i32) -> i32 {
        self.seed = (self.seed.wrapping_mul(0x5DEECE66D).wrapping_add(0xB))
            & ((1u64 << 48) - 1);
        (self.seed >> (48 - bits)) as i32
    }

    pub fn next_int(&mut self, n: i32) -> i32 {
        if n <= 0 {
            panic!("n must be positive");
        }
        if (n & -n) == n {
            return ((n as u64 * self.next(31) as u64) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let val = bits % n;
            if bits - val + (n - 1) >= 0 {
                return val;
            }
        }
    }
}

pub fn generate_random_permutation(n_elements: i32, random: &mut JavaRandom) -> Vec<i32> {
    let mut permutation: Vec<i32> = (0..n_elements).collect();
    for i in 0..n_elements {
        let j = random.next_int(n_elements);
        let k = permutation[i as usize];
        permutation[i as usize] = permutation[j as usize];
        permutation[j as usize] = k;
    }
    permutation
}
