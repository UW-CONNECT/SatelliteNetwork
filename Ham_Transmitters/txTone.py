#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# GNU Radio version: 3.10.9.2

from PyQt5 import Qt
from gnuradio import qtgui
from PyQt5 import QtCore
from gnuradio import analog
from gnuradio import gr
from gnuradio.filter import firdes
from gnuradio.fft import window
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time
import sip



class txTone(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Not titled yet", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Not titled yet")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except BaseException as exc:
            print(f"Qt GUI: Could not set Icon: {str(exc)}", file=sys.stderr)
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "txTone")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.tx_power_dB = tx_power_dB = 10
        self.samp_rate = samp_rate = 192e3

        ##################################################
        # Blocks
        ##################################################

        self._tx_power_dB_range = qtgui.Range(0, 89, 10, 10, 200)
        self._tx_power_dB_win = qtgui.RangeWidget(self._tx_power_dB_range, self.set_tx_power_dB, "'tx_power_dB'", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_layout.addWidget(self._tx_power_dB_win)
        self.uhd_usrp_sink_0_0 = uhd.usrp_sink(
            ",".join(("", "")),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
            '',
        )
        self.uhd_usrp_sink_0_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0_0.set_time_unknown_pps(uhd.time_spec(0))

        self.uhd_usrp_sink_0_0.set_center_freq(435e6, 0)
        self.uhd_usrp_sink_0_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_sink_0_0.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_sink_0_0.set_gain(tx_power_dB, 0)
        self.qtgui_sink_x_0_0 = qtgui.sink_c(
            1024, #fftsize
            window.WIN_BLACKMAN_hARRIS, #wintype
            0, #fc
            samp_rate, #bw
            "", #name
            True, #plotfreq
            True, #plotwaterfall
            True, #plottime
            True, #plotconst
            None # parent
        )
        self.qtgui_sink_x_0_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_0_win = sip.wrapinstance(self.qtgui_sink_x_0_0.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_0_0.enable_rf_freq(False)

        self.top_layout.addWidget(self._qtgui_sink_x_0_0_win)
        self.analog_sig_source_x_0_0 = analog.sig_source_c(samp_rate, analog.GR_SIN_WAVE, (-1e3), 1, 0, 0)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_sig_source_x_0_0, 0), (self.qtgui_sink_x_0_0, 0))
        self.connect((self.analog_sig_source_x_0_0, 0), (self.uhd_usrp_sink_0_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "txTone")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_tx_power_dB(self):
        return self.tx_power_dB

    def set_tx_power_dB(self, tx_power_dB):
        self.tx_power_dB = tx_power_dB
        self.uhd_usrp_sink_0_0.set_gain(self.tx_power_dB, 0)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.analog_sig_source_x_0_0.set_sampling_freq(self.samp_rate)
        self.qtgui_sink_x_0_0.set_frequency_range(0, self.samp_rate)
        self.uhd_usrp_sink_0_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_0_0.set_bandwidth(self.samp_rate, 0)




def main(top_block_cls=txTone, options=None):

    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    qapp.exec_()

if __name__ == '__main__':
    main()
