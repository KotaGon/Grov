COMPILER = clang++
CFLAGS   = -g -Wall -O3 -std=c++17
LDFLAGS  = 
LIBS     = 
INCLUDE  = -I./include -I./
TARGET   = optimize
OBJDIR   = ./obj
INCDIR   = ./include
SOURCES  = $(wildcard *.cc *.cpp)
HEADERS  = $(wildcard $(INCDIR)/*.h $(INCDIR)/*.hpp)
OBJECTS  = $(patsubst %.cc, $(OBJDIR)/%.o, $(filter %.cc, $(SOURCES))) \
           $(patsubst %.cpp, $(OBJDIR)/%.o, $(filter %.cpp, $(SOURCES)))

$(TARGET): $(OBJECTS) $(LIBS)
	$(COMPILER) -o $@ $^ $(LDFLAGS)
	
$(OBJDIR)/%.o: %.cc
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(OBJDIR)/%.o: %.cpp
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(OBJECTS): $(HEADERS)

all: clean $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)
