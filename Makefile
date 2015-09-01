SRC_DIR = ./src
OBJ_DIR = ./obj
lapack_SRC_DIR = ./src_Lapack
lapack_OBJ_DIR = ./obj_Lapack
PROG = ForwardIBD 
OBJS = $(OBJ_DIR)/random.o $(OBJ_DIR)/pops.o $(OBJ_DIR)/params.o $(OBJ_DIR)/matrix.o $(OBJ_DIR)/expansion.o $(lapack_OBJ_DIR)/*.o 
CC = gcc
CFLAGS = -g -lm -O3 -Wformat=0 -Wno-div-by-zero 

$(PROG): $(OBJS) $(OBJS_c)
	$(CC) $(SRC_DIR)/Forward.c -o $@ $(OBJS) $(CFLAGS)

$(OBJ_DIR)/random.o: $(SRC_DIR)/random.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/random.c -o $(OBJ_DIR)/random.o

$(OBJ_DIR)/pops.o: $(SRC_DIR)/pops.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/pops.c -o $(OBJ_DIR)/pops.o

$(OBJ_DIR)/params.o: $(SRC_DIR)/params.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/params.c -o $(OBJ_DIR)/params.o

$(OBJ_DIR)/matrix.o: $(SRC_DIR)/matrix.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/matrix.c -o $(OBJ_DIR)/matrix.o

$(OBJ_DIR)/expansion.o: $(SRC_DIR)/expansion.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/expansion.c -o $(OBJ_DIR)/expansion.o

lapack: 
	$(CC) $(CFLAGS) -c $(lapack_SRC_DIR)/*.c
	mv *.o $(lapack_OBJ_DIR)/

clean:
	rm -f $(OBJ_DIR)/*.o $(PROG) 

realclean:
	rm -f $(OBJ_DIR)/*.o $(PROG) 
	rm -f $(OBJ__c_DIR)/*.o $(PROG_c1) $(PROG_c2)
	rm -f $(lapack_OBJ_DIR)/*.o $(OBJ_DIR)/*.o $(PROG) 

