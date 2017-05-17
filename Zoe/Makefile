#############################
# Makefile for ZOE library  #
#############################

LIB = -lm

OBJECTS = \
	zoeCDS.o\
	zoeCounter.o\
	zoeDNA.o\
	zoeDistribution.o\
	zoeDuration.o\
	zoeFastaFile.o\
	zoeFeature.o\
	zoeFeatureFactory.o\
	zoeFeatureTable.o\
	zoeHMM.o\
	zoeIsochore.o\
	zoeMath.o\
	zoeModel.o\
	zoePhasePref.o\
	zoeProtein.o\
	zoeScanner.o\
	zoeState.o\
	zoeTools.o\
	zoeTransition.o\
	zoeTrellis.o\


APP = zoe-loop
SRC = zoe-loop.c
OBJ = zoe-loop.o

DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(CFLAGS) $(OBJ) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP) $(APP2) $(APP3) $(APP4) $(APP5)

depend: $(OBJECTS:.o=.c)
	gcc -MM $^ > $@

tar:
	rm -rf /tmp/$(APP)
	mkdir /tmp/$(APP)
	cp -r * /tmp/$(APP)
	cd /tmp/$(APP); make clean;  rm -rf CVS */CVS
	cd /tmp; tar -zcvf $(APP)-$(DATE).tar.gz $(APP)
	rm -rf /tmp/$(APP)


#################
# Architectures #
#################

gcc:
	make $(APP)  CC="gcc" CFLAGS="-O2 -Wall -Werror"

###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

################
# Dependancies #
################

include depend

