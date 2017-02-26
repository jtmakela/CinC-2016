BASE := .

SRCDIR := ./src
OBJDIR := $(BASE)/obj
BINDIR := $(BASE)/bin

CC := g++
LD := g++

CC_FLAGS := -Wall -march=native -O3 -I "./include"
LD_FLAGS :=
LD_LIBS := -pthread -lm

all: $(BASE) $(BINDIR) $(OBJDIR) $(BINDIR)/tftrig_final 
	
$(BINDIR)/tftrig_final: $(SRCDIR)/tftrig_final.cpp $(OBJDIR)/Simplified.so $(OBJDIR)/myDSP.so $(OBJDIR)/Trigger.so $(OBJDIR)/Classifier.so
	$(CC) $(CC_FLAGS) -o $@ $^ $(LD_LIBS)

$(BASE):
	if [ ! -d $(BASE) ]; then \
		sudo mkdir -m 0755 $(BASE);\
		sudo chown $(shell whoami): $(BASE);\
	fi

$(OBJDIR):
	if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/Simplified.so: $(SRCDIR)/Simplified/PhysionetChallenge2016.cpp $(SRCDIR)/Simplified/Retrigger.cpp
	$(CC) $(CC_FLAGS) -o $@ -shared -fPIC $^

$(OBJDIR)/myDSP.so: $(SRCDIR)/myDSP/Iir.cpp
	$(CC) $(CC_FLAGS) -o $@ -shared -fPIC $^

$(OBJDIR)/Classifier.so: $(SRCDIR)/Classifier/classifier.cpp $(SRCDIR)/Classifier/markers.cpp
	$(CC) $(CC_FLAGS) -o $@ -shared -fPIC $^

$(OBJDIR)/Trigger.so: $(SRCDIR)/Trigger/Csv2kernel.cpp $(SRCDIR)/Trigger/Trigger.cpp
	$(CC) $(CC_FLAGS) -o $@ -shared -fPIC $^

$(BINDIR):
	if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi

clean:
	rm $(BINDIR)/tftrig_final

entry:	
	$(shell find . -type d -name CVS | xargs rm -r)
	if [ -e entry.zip ]; then rm entry.zip; fi
	zip -r entry.zip AUTHORS.txt LICENSE.txt src include Makefile next.sh setup.sh trigger.convolution.kernel.csv answers.txt 

dryrun:	
	$(shell find . -type d -name CVS | xargs rm -r)
	if [ -e entry.zip ]; then rm entry.zip; fi
	zip -r entry.zip DRYRUN AUTHORS.txt LICENSE.txt src include Makefile next.sh setup.sh trigger.convolution.kernel.csv answers.txt 
