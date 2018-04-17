CC          = icpc
CLINKER     = icpc

#CFLAGS      =   -Wall -O3 -march=pentium
CFLAGS      =   -Wall -O3 -xHost

#CFLAGS      = -i-fast -lm  
#LIBS        = -lm
DEPEND= makedepend

SRC        = PatchyMD.c ran_uniform.c system.c boxmuller.c randomsphere.c readinput.c initialization.c write.c potential.c nve.c nve2.c nvt.c transfer.c force.c kinetic.c gr.c scale.c 
OBJS       = PatchyMD.o ran_uniform.o system.o boxmuller.o randomsphere.o readinput.o initialization.o write.o potential.o nve.o nve2.c nvt.o transfer.o force.o kinetic.o gr.o scale.o
EXECS      = PatchyMD

default: PatchyMD

all: $(EXECS)

PatchyMD:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o PatchyMD $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

PatchyMD.o: system.h ran_uniform.h 
ran_uniform.o: system.h ran_uniform.h
boxmuller.o: system.h ran_uniform.h
randomsphere.o: system.h ran_uniform.h
system.o: system.h
readinput.o: system.h ran_uniform.h
initialization.o: system.h ran_uniform.h
write.o: system.h
potential.o: system.h
nve.o: system.h
nve2.o: system.h
nvt.o: system.h
transfer.o: system.h
force.o: system.h
kinetic.o: system.h
gr.o: system.h
scale.o: system.h
