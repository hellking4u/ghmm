# ghmm/doc/Makefile.am

#******************  Doc++ *******************

DOCXX            = @DOCXX@
DVIPS            = @DVIPS@
LATEX            = @LATEX@

#Name der latex Dokumentation
APP              = @PACKAGE@
#Verzeichnisse für HTML und latex Dokumentation
HTML             = docHTML
DOC              = docTEX

#eine sammelnde Doku-Datei mit Includes auf die .h-Files
DOCXX_FILES      = docu.dxx
BUILT_SOURCES   += $(DOCXX_FILES)

# html installation
htmldir   = $(prefix)/html
html_DATA = docHTML

install-htmlDATA:
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(htmldir)/$(APP)
	@list='docHTML/*'; for p in $$list; do \
	  if test -f $$p; then \
	    f=`basename $$p` ;\
	    echo " $(INSTALL_DATA) $$p $(DESTDIR)$(htmldir)/hmmlib/$$f"; \
	    $(INSTALL_DATA) $$p $(DESTDIR)$(htmldir)/$(APP)/$$f; \
	  else :; fi; \
	done
	@$(POST_INSTALL)

uninstall-htmlDATA:
	@$(NORMAL_UNINSTALL)
	list='$(html_DATA)'; for p in $$list; do \
	  f=`basename $$p` ;\
	  rm -f $(DESTDIR)$(htmldir)/$(APP)/$$f; \
	done
	-rmdir $(DESTDIR)$(htmldir)/$(APP)

# Wahrscheinlich nicht zu editierende Optionen:

DOCXX_OPTS       = -a -p
DOCXX_LATEX_OPTS = -t -eo a4paper
DOCXX_HTML_OPTS  = #-H

#for make (dist)clean
CLEANFILES += docHTML/* docTEX/$(APP).*
EXTRA_DIST += docTEX/docxx.sty

$(HTML)/index.html: $(DOCXX_FILES)
	test -z "$(DOCXX)" || $(DOCXX) $(DOCXX_OPTS) $(DOCXX_HTML_OPTS) -d $(HTML) $(DOCXX_FILES)

$(DOC)/$(APP).tex: $(DOCXX_FILES)
	test -z "$(DOCXX)" || $(DOCXX) $(DOCXX_OPTS) $(DOCXX_LATEX_OPTS) -o $(DOC)/$(APP).tex $(DOCXX_FILES)

$(DOC)/$(APP).dvi: $(DOC)/$(APP).tex
	cd $(DOC); $(LATEX) $(APP).tex; latex $(APP).tex   

$(DOC)/$(APP).ps: $(DOC)/$(APP).dvi
	test -z "$(DVIPS)" || $(DVIPS) $(DOC)/$(APP).dvi -o $(DOC)/$(APP).ps

docTEX: $(DOC)/$(APP).ps

docHTML: $(HTML)/index.html

