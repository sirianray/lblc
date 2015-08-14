CC = icc
CFLAG = -O3
SOURCE = main.c lb.h func.c init.c lattice.c read_param.c particle.h particle.c fd_Q.c evolution.c
TARGET = lblc

$(TARGET): $(SOURCE)
	$(CC) $(CFLAG) -o $(TARGET) $(SOURCE)

clean:
	rm $(TARGET) 
