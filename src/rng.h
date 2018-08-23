

/* Convert u64 to uniform f64. See http://xoshiro.di.unimi.it.
This will cut the number of possible values in half as the lowest significant
bit will be set to 0 for all returned values. */

Inline
f64 toFloat(u64 x)
{
    const union { u64 i; f64 d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
}


// NOTE(jonas): rng abstraction might be too slow because inlining not possible

#define NEXT_FUNC(name) Inline u64 name(void *state)
typedef u64 nextFunc(void *state);

#define NEXT_FLOAT_FUNC(name) Inline f64 name(void *state)
typedef f64 nextFloatFunc(void *state);

#define JUMP_FUNC(name) void name(void *state)
typedef void jumpFunc(void *state);


typedef struct Rng Rng;
struct Rng {
    nextFunc      *next;
    nextFloatFunc *nextFloat;
    jumpFunc      *jump;
    void *state;
};


Inline
u64 RngNext( Rng r ) {
    return r.next( r.state );
}

Inline
f64 RngNextFloat( Rng r ) {
    return r.nextFloat( r.state );
}

Inline
void RngJump( Rng r ) {
    r.jump( r.state );
}



/* xorshift1024*, see http://vigna.di.unimi.it/ftp/papers/xorshift.pdf */

typedef struct Xorshift1024 Xorshift1024;
struct Xorshift1024 {
    u64 s[16];
    i32 p;
};


Xorshift1024 Xorshift1024Init( u64 seed )
{
    Xorshift1024 state;

    for ( u32 i=0; i<16; ++i ) {
        state.s[i] = seed + i;
    }

    state.p = 0;

    return state;
}

#define XORSHIFT_INIT Xorshift1024Init( (u64) time(NULL) )


NEXT_FUNC(rngXorshift1024Next) {

    Xorshift1024 *x = (Xorshift1024 *) state;

    const u64 s0 = x->s[x->p];
    u64       s1 = x->s[ x->p = (x->p + 1) & 15 ];

    s1 ^= s1 << 31;

    x->s[x->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);

    return x->s[x->p] * UINT64_C(1181783497276652981);
}


NEXT_FLOAT_FUNC(rngXorshift1024NextFloat)
{
    u64 x = rngXorshift1024Next( state );

    return toFloat( x );
}


JUMP_FUNC(rngXorshift1024Jump)
{
    Xorshift1024 *x = (Xorshift1024 *) state;

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
                    t[j] ^= x->s[(j + x->p) & 15];
                }
            }
            rngXorshift1024Next( state );
        }
    }

    for ( i32 j = 0; j < 16; j++ ) {
        x->s[(j + x->p) & 15] = t[j];
    }
}


Rng RngInitXorshift1024( Xorshift1024 *state )
{
    Rng r;
    r.next      = rngXorshift1024Next;
    r.nextFloat = rngXorshift1024NextFloat;
    r.jump      = rngXorshift1024Jump;
    r.state     = state;
    return r;
}


#if TEST
void test_xorshift1024()
{
    Xorshift1024 x = Xorshift1024Init( 37473 );
    
    TEST_ASSERT( rngXorshift1024Next( &x ) > 0 );
    TEST_ASSERT( rngXorshift1024NextFloat( &x ) < 1.0 );
    
    rngXorshift1024Jump( &x);
    
    TEST_ASSERT( rngXorshift1024Next( &x ) > 0 );
    TEST_ASSERT( rngXorshift1024NextFloat( &x ) < 1.0 );
}
#endif


// NOTE(jonas): macros to allow inlining

#define RNG_TYPE Xorshift1024
#define RNG_NEXT rngXorshift1024Next
#define RNG_NEXT_FLOAT rngXorshift1024NextFloat
#define RNG_JUMP rngXorshift1024Jump



/* Box-Muller Transformation for normal distribution samples.
 See: https://en.wikipedia.org/wiki/Box–Muller_transform .*/

Inline
f64 RngNormal( Rng r ){

    f64 u = RngNextFloat( r );
    f64 v = RngNextFloat( r );

    f64 s = sqrt( -2 * log( u ) );
    f64 t = cos( 2 * M_PI * v );


    return s * t;
}



// TODO(jonas): use exact normal sampling

/* Exact sampling from normal distribution: https://arxiv.org/abs/1303.6257 .
 The algorithm is exact when the underlying uniform RNG is perfect. */

