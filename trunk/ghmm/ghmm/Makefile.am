# ghmm/lib/Makefile.am

lib_LIBRARIES = libghmm.a

libghmm_a_SOURCES= foba.c foba.h \
                  matrix.c  matrix.h \
                  model.c model.h \
                  cluster.c cluster.h \
                  scluster.c scluster.h \
		  root_finder.c root_finder.h \
                  reestimate.c reestimate.h \
                  randvar.c randvar.h \
                  rng.c rng.h \
                  sequence.c sequence.h \
                  vector.c vector.h \
                  viterbi.c viterbi.h \
                  sfoba.c sfoba.h \
                  sgenerate.c sgenerate.h \
                  smodel.c smodel.h \
                  sreestimate.c sreestimate.h \
                  smixturehmm.c smixturehmm.h \
                  const.h \
                  mes.c mes.h \
                  mprintf.c mprintf.h \
                  scanner.c scanner.h \
                  gauss_tail.c gauss_tail.h


include_HEADERS = foba.h \
                  matrix.h \
                  model.h \
                  cluster.h \
                  scluster.h \
		  root_finder.h \
                  reestimate.h \
                  randvar.h \
                  rng.h \
                  sequence.h \
                  vector.h \
                  viterbi.h \
                  sfoba.h \
                  sgenerate.h \
                  smodel.h \
                  sreestimate.h \
                  smixturehmm.h \
                  const.h \
                  mes.h \
                  mprintf.h \
                  scanner.h \
                  gauss_tail.h

install-includeHEADERS: $(include_HEADERS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(includedir)
	$(mkinstalldirs) $(DESTDIR)$(includedir)/ghmm
	@list='$(include_HEADERS)'; for p in $$list; do \
	  if test -f "$$p"; then d= ; else d="$(srcdir)/"; fi; \
	  f=`basename $$p`; \
	  echo " $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/ghmm/$$f"; \
	 $(INSTALL_DATA) $$d$$p $(DESTDIR)$(includedir)/ghmm/$$f; \
	done

uninstall-includeHEADERS:
	@$(NORMAL_UNINSTALL)
	list='$(include_HEADERS)'; for p in $$list; do \
	  f=`basename $$p`; \
	  rm -f $(DESTDIR)$(includedir)/ghmm/$$f; \
	done
	-rmdir $(DESTDIR)$(includedir)/ghmm
