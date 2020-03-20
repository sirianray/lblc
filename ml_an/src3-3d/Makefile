CC = mpicc
CFLAG = -O3 -mkl
SOURCE = main.c main.h read_param.c func.c init.c memory.c geometry.c fd_Q.c evolution.c particle.h particle.c patches.c
TARGET = lblc

$(TARGET): $(SOURCE)
	$(CC) $(CFLAG) -o $(TARGET) $(SOURCE)

clean:
	rm $(TARGET) *.out 
