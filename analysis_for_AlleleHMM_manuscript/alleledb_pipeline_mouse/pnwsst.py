#
# Copyright (c) 2005-2008, REvolution Computing, Inc.
#
# All rights reserved
#

import socket, sys, threading
from twisted.application import internet
from twisted.internet import reactor
from twisted.python import log
from nwss.server import NwsService
try:
    from twisted.web import server
except ImportError:
    server = None

__all__ = ['NwsLocalServerException', 'NwsLocalServer']

class NwsLocalServerException(Exception):
    pass

class NwsLocalServer(threading.Thread):
    invoked = False

    def __init__(self, port=0, interface='', daemon=True,
                 name='NwsLocalServer', logFile=None, quiet=False, **kw):

        if logFile:
            log.startLogging(open(logFile, 'w'), setStdout=False)
        elif not quiet:
            log.startLogging(sys.stderr, setStdout=False)

        if NwsLocalServer.invoked:
            raise NwsLocalServerException('I can only be invoked once!')
        NwsLocalServer.invoked = True

        threading.Thread.__init__(self, name=name, **kw)

        self._desiredPort = port
        self._interface = interface
        self.setDaemon(daemon)
        self.startSem = threading.Semaphore(0)
        self.stopSem = threading.Semaphore(0)
        self.stoppedSem = threading.Semaphore(0)
        self.start()
        self.startSem.acquire()

    def run(self):
        srv = NwsService(self._desiredPort, interface=self._interface)
        srv.startService()  # XXX do this after listenTCP?
        self._p = reactor.listenTCP(self._desiredPort, srv._factory)
        self.host = socket.gethostname() # is there some field in self._p that has this info?
        self.port = self._p._realPortNumber
        if server:
            self._wp = reactor.listenTCP(0, server.Site(srv.getResource()))
            self.webport = self._wp._realPortNumber
        else:
            self._wp = None
            self.webport = -1
        self.startSem.release()
        def shutDown():
            self.stopSem.acquire()
            log.msg('Got stop semaphore. Adding call to reactor.stop.')
            reactor.callFromThread(reactor.stop)
        reactor.callInThread(shutDown)
        try:
            reactor.run(installSignalHandlers=0)
        finally:
            srv.stopService()
            self.stoppedSem.release()

    def stop(self):
        self.stopSem.release()
        self.stoppedSem.acquire()
