

/* Convert u64 to uniform f64. See http://xoshiro.di.unimi.it.
This will cut the number of possible values in half as the lowest significant
bit will be set to 0 for all returned values. */

Inline
f64 toFloat(u64 x)
{
    const union { u64 i; f64 d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
}


/* xorshift1024*, see http://vigna.di.unimi.it/ftp/papers/xorshift.pdf */

typedef struct StateXorshift1024 StateXorshift1024;
struct StateXorshift1024 {
    u64 s[16];
    i32 p;
};


u64 xorshift1024( StateXorshift1024 *state )
{
    const u64 s0 = state->s[state->p];
    u64       s1 = state->s[ state->p = (state->p + 1) & 15 ];
    
    s1 ^= s1 << 31;
    
    state->s[state->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);
    
    return state->s[state->p] * UINT64_C(1181783497276652981);
}


void jump( StateXorshift1024 *state )
{
    static const u64 JUMP[] = {
        0x84242f96eca9c41d, 0xa3c65b8776f96855, 0x5b34a39f070b5837,
        0x4489affce4f31a1e, 0x2ffeeb0a48316f40, 0xdc2d9891fe68c022,
        0x3659132bb12fea70, 0xaac17d8efa43cab8, 0xc4cb815590989b13,
        0x5ee975283d71c93b, 0x691548c86c1bd540, 0x7910c41d10a1e6a5,
        0x0b5fc64563b3e2a8, 0x047f7684e9fc949d, 0xb99181f2d8f685ca,
        0x284600e3f30e38c3
    };
    
    u64 t[16] = { 0 };
    
    for ( i32 i = 0; i < ( sizeof(JUMP) / sizeof(*JUMP) ); ++i ) {
        for (i32 b = 0; b < 64; ++b) {
            if ( JUMP[i] & 1ULL << b ) {
                for (i32 j = 0; j < 16; ++j ) {
                    t[j] ^= state->s[(j + state->p) & 15];
                }
            }
            xorshift1024( state );
        }
    }
    
    for ( i32 j = 0; j < 16; j++ ) {
        state->s[(j + state->p) & 15] = t[j];
    }
}


f64 floatXorshift1024( StateXorshift1024 *state )
{
    u64 x = xorshift1024( state );
    
    return toFloat( x );
}

