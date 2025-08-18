#ifndef RANDOM_H
#define RANDOM_H

#define RAND48_SEED_0 (0x330e)
#define RAND48_SEED_1 (0xabcd)
#define RAND48_SEED_2 (0x1234)
#define RAND48_MULT_0 (0xe66d)
#define RAND48_MULT_1 (0xdeec)
#define RAND48_MULT_2 (0x0005)
#define RAND48_ADD (0x000b)

struct random_state {
    random_state(int seed) {
        x[0] = RAND48_SEED_0;
        x[1] = (unsigned short)seed;
        x[2] = (unsigned short)(seed >> 16);
    };

    unsigned short x[3];
};

extern double my_erand48(unsigned short xseed[3]);
extern double my_drand48();
extern void my_srand48(long seed);

#endif