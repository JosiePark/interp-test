include ../../standard_defs.mk

all: *.o move


*.o : *.f
	$(CC) -c $(LIBS) $(INCLUDES) *.f $(LINKS)
move:
	mv *.o *.mod $(build_DIR)

clean:
	rm *.o *.mod 	

