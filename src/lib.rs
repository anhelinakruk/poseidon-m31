const P: u32 = (1u32 << 31) - 1;
const N_STATE: usize = 16;
const RATE: usize = 8;
const N_HALF_FULL_ROUNDS: usize = 4;
const N_PARTIAL_ROUNDS: usize = 14;

const MAT_INTERNAL_DIAG_M_1: [u32; N_STATE] = [
    0x07b80ac4, 0x6bd9cb33, 0x48ee3f9f, 0x4f63dd19,
    0x18c546b3, 0x5af89e8b, 0x4ff23de8, 0x4f78aaf6,
    0x53bdc6d4, 0x5c59823e, 0x2a471c72, 0x4c975e79,
    0x58dc64d4, 0x06e9315d, 0x2cf32286, 0x2fb6755d,
];

const EXTERNAL_ROUND_CONSTS: [[u32; N_STATE]; 2 * N_HALF_FULL_ROUNDS] = [
    [0x768bab52, 0x70e0ab7d, 0x3d266c8a, 0x6da42045, 0x600fef22, 0x41dace6b, 0x64f9bdd4, 0x5d42d4fe, 0x76b1516d, 0x6fc9a717, 0x70ac4fb6, 0x00194ef6, 0x22b644e2, 0x1f7916d5, 0x47581be2, 0x2710a123],
    [0x6284e867, 0x018d3afe, 0x5df99ef3, 0x4c1e467b, 0x566f6abc, 0x2994e427, 0x538a6d42, 0x5d7bf2cf, 0x7fda2dab, 0x0fd854c4, 0x46922fca, 0x3d7763a1, 0x19fd05ca, 0x0a4bbb43, 0x15075851, 0x3d903d76],
    [0x2d290ff7, 0x40809fa0, 0x59dac6ec, 0x127927a2, 0x6bbf0ea0, 0x0294140f, 0x24742976, 0x6e84c081, 0x22484f4a, 0x354cae59, 0x0453ffe1, 0x3f47a3cc, 0x0088204e, 0x6066e109, 0x3b7c4b80, 0x6b55665d],
    [0x3bc4b897, 0x735bf378, 0x508daf42, 0x1884fc2b, 0x7214f24c, 0x7498be0a, 0x1a60e640, 0x3303f928, 0x29b46376, 0x5c96bb68, 0x65d097a5, 0x1d358e9f, 0x4a9a9017, 0x4724cf76, 0x347af70f, 0x1e77e59a],
    [0x57090613, 0x1fa42108, 0x17bbef50, 0x1ff7e11c, 0x047b24ca, 0x4e140275, 0x4fa086f5, 0x079b309c, 0x1159bd47, 0x6d37e4e5, 0x075d8dce, 0x12121ca0, 0x7f6a7c40, 0x68e182ba, 0x5493201b, 0x0444a80e],
    [0x0064f4c6, 0x6467abe6, 0x66975762, 0x2af68f9b, 0x345b33be, 0x1b70d47f, 0x053db717, 0x381189cb, 0x43b915f8, 0x20df3694, 0x0f459d26, 0x77a0e97b, 0x2f73e739, 0x1876c2f9, 0x65a0e29a, 0x4cabefbe],
    [0x5abd1268, 0x4d34a760, 0x12771799, 0x69a0c9ac, 0x39091e55, 0x7f611cd0, 0x3af055da, 0x7ac0bbdf, 0x6e0f3a24, 0x41e3b6f7, 0x49b3756d, 0x568bc538, 0x20c079d8, 0x1701c72c, 0x7670dc6c, 0x5a439035],
    [0x7c93e00e, 0x561fbb4d, 0x1178907b, 0x02737406, 0x32fb24f1, 0x6323b60a, 0x6ab12418, 0x42c99cea, 0x155a0b97, 0x53d1c6aa, 0x2bd20347, 0x279b3d73, 0x4f5f3c70, 0x0245af6c, 0x238359d3, 0x49966a59],
];

const INTERNAL_ROUND_CONSTS: [u32; N_PARTIAL_ROUNDS] = [
    0x7f7ec4bf, 0x0421926f, 0x5198e669, 0x34db3148, 0x4368bafd, 0x66685c7f,
    0x78d3249a, 0x60187881, 0x76dad67a, 0x0690b437, 0x1ea95311, 0x40e5369a,
    0x38f103fc, 0x1d226a21,
];

#[inline]
fn add(a: u32, b: u32) -> u32 {
    let s = a.wrapping_add(b);
    if s >= P { s - P } else { s }
}

#[inline]
fn mul(a: u32, b: u32) -> u32 {
    ((a as u64 * b as u64) % P as u64) as u32
}

#[inline]
fn pow5(x: u32) -> u32 {
    let x2 = mul(x, x);
    let x4 = mul(x2, x2);
    mul(x4, x)
}


fn apply_m4(x: [u32; 4]) -> [u32; 4] {
    let t0 = add(x[0], x[1]);
    let t1 = add(x[2], x[3]);
    let t02 = add(t0, t0);
    let t12 = add(t1, t1);
    let x1d = add(x[1], x[1]);
    let t2 = add(x1d, t1);
    let x3d = add(x[3], x[3]);
    let t3 = add(x3d, t0);
    let t4 = add(add(t12, t12), t3);
    let t5 = add(add(t02, t02), t2);
    [add(t3, t5), t5, add(t2, t4), t4]
}

fn external_matrix(mut s: [u32; N_STATE]) -> [u32; N_STATE] {
    for i in 0..4 {
        let c = apply_m4([s[4*i], s[4*i+1], s[4*i+2], s[4*i+3]]);
        s[4*i..4*i+4].copy_from_slice(&c);
    }
    for j in 0..4 {
        let t = add(add(s[j], s[j+4]), add(s[j+8], s[j+12]));
        s[j]    = add(s[j],    t);
        s[j+4]  = add(s[j+4],  t);
        s[j+8]  = add(s[j+8],  t);
        s[j+12] = add(s[j+12], t);
    }
    s
}

fn internal_matrix(mut s: [u32; N_STATE]) -> [u32; N_STATE] {
    let sum = s.iter().fold(0u32, |acc, &x| add(acc, x));
    for i in 0..N_STATE {
        s[i] = add(mul(s[i], MAT_INTERNAL_DIAG_M_1[i]), sum);
    }
    s
}

fn permute(mut s: [u32; N_STATE]) -> [u32; N_STATE] {
    s = external_matrix(s);
    for round in 0..N_HALF_FULL_ROUNDS {
        for i in 0..N_STATE { s[i] = add(s[i], EXTERNAL_ROUND_CONSTS[round][i]); }
        for i in 0..N_STATE { s[i] = pow5(s[i]); }
        s = external_matrix(s);
    }
    for round in 0..N_PARTIAL_ROUNDS {
        s[0] = pow5(add(s[0], INTERNAL_ROUND_CONSTS[round]));
        s = internal_matrix(s);
    }
    for round in N_HALF_FULL_ROUNDS..2 * N_HALF_FULL_ROUNDS {
        for i in 0..N_STATE { s[i] = add(s[i], EXTERNAL_ROUND_CONSTS[round][i]); }
        for i in 0..N_STATE { s[i] = pow5(s[i]); }
        s = external_matrix(s);
    }
    s
}

fn absorb(mut s: [u32; N_STATE], block: [u32; RATE]) -> [u32; N_STATE] {
    for i in 0..RATE { s[i] = add(s[i], block[i]); }
    permute(s)
}

/// Streaming Poseidon2 hasher over M31.
///
/// # Example
/// ```
/// use posiedon_m31::PoseidonHasher;
///
/// let mut h = PoseidonHasher::new();
/// h.update(b"hello ");
/// h.update(b"world");
/// let digest = h.finalize();
/// ```
pub struct PoseidonHasher {
    state: [u32; N_STATE],
    buf: [u32; RATE],
    pos: usize,
    absorbed: bool,
}

impl PoseidonHasher {
    pub fn new() -> Self {
        Self {
            state: [0u32; N_STATE],
            buf: [0u32; RATE],
            pos: 0,
            absorbed: false,
        }
    }

    /// Feed more bytes into the hasher.
    pub fn update(&mut self, data: &[u8]) {
        for &b in data {
            self.buf[self.pos] = b as u32;
            self.pos += 1;
            if self.pos == RATE {
                self.state = absorb(self.state, self.buf);
                self.buf = [0u32; RATE];
                self.pos = 0;
                self.absorbed = true;
            }
        }
    }

    /// Consume the hasher and return the 32-byte digest.
    pub fn finalize(mut self) -> [u8; 32] {
        if !self.absorbed || self.pos > 0 {
            self.state = absorb(self.state, self.buf);
        }
        let mut out = [0u8; 32];
        for (i, &w) in self.state[..RATE].iter().enumerate() {
            out[i * 4..(i + 1) * 4].copy_from_slice(&w.to_le_bytes());
        }
        out
    }
}

impl Default for PoseidonHasher {
    fn default() -> Self {
        Self::new()
    }
}

/// Hashes a byte slice with Poseidon2 over M31.
pub fn hash_bytes(data: &[u8]) -> [u8; 32] {
    let mut h = PoseidonHasher::new();
    h.update(data);
    h.finalize()
}
