#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: NFM_xmt
# Author: Barry Duggan
# Description: NBFM transmitter
# GNU Radio version: 3.10.9.2

from PyQt5 import Qt
from gnuradio import qtgui
from PyQt5 import QtCore
from PyQt5.QtCore import QObject, pyqtSlot
from gnuradio import analog
from gnuradio import blocks
import pmt
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
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



class NFM_xmt(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "NFM_xmt", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("NFM_xmt")
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

        self.settings = Qt.QSettings("GNU Radio", "NFM_xmt")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.usrp_rate = usrp_rate = 576000
        self.volume_gain = volume_gain = 1
        self.samp_rate = samp_rate = 48000
        self.ptt = ptt = 1
        self.pl_freq = pl_freq = 67
        self.if_rate = if_rate = int(usrp_rate/3)
        self.freq_tx = freq_tx = 145.99

        ##################################################
        # Blocks
        ##################################################

        self._volume_gain_range = qtgui.Range(0, 100, 0.1, 1, 200)
        self._volume_gain_win = qtgui.RangeWidget(self._volume_gain_range, self.set_volume_gain, "Audio gain", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_layout.addWidget(self._volume_gain_win)
        _ptt_push_button = Qt.QPushButton('Push To Talk')
        _ptt_push_button = Qt.QPushButton('Push To Talk')
        self._ptt_choices = {'Pressed': 1, 'Released': 0}
        _ptt_push_button.pressed.connect(lambda: self.set_ptt(self._ptt_choices['Pressed']))
        _ptt_push_button.released.connect(lambda: self.set_ptt(self._ptt_choices['Released']))
        self.top_layout.addWidget(_ptt_push_button)
        # Create the options list
        self._pl_freq_options = [0.0, 67.0, 71.9, 74.4, 77.0, 79.7, 82.5, 85.4, 88.5, 91.5, 94.8, 97.4, 100.0, 103.5, 107.2, 110.9, 114.8, 118.8, 123.0, 127.3, 131.8, 136.5, 141.3, 146.2, 151.4, 156.7, 162.2, 167.9, 173.8, 179.9, 186.2, 192.8, 203.5, 210.7, 218.1, 225.7, 233.6, 241.8, 250.3]
        # Create the labels list
        self._pl_freq_labels = ['0.0', '67.0', '71.9', '74.4', '77.0', '79.7', '82.5', '85.4', '88.5', '91.5', '94.8', '97.4', '100.0', '103.5', '107.2', '110.9', '114.8', '118.8', '123.0', '127.3', '131.8', '136.5', '141.3', '146.2', '151.4', '156.7', '162.2', '167.9', '173.8', '179.9', '186.2', '192.8', '203.5', '210.7', '218.1', '225.7', '233.6', '241.8', '250.3']
        # Create the combo box
        self._pl_freq_tool_bar = Qt.QToolBar(self)
        self._pl_freq_tool_bar.addWidget(Qt.QLabel("PL Tone" + ": "))
        self._pl_freq_combo_box = Qt.QComboBox()
        self._pl_freq_tool_bar.addWidget(self._pl_freq_combo_box)
        for _label in self._pl_freq_labels: self._pl_freq_combo_box.addItem(_label)
        self._pl_freq_callback = lambda i: Qt.QMetaObject.invokeMethod(self._pl_freq_combo_box, "setCurrentIndex", Qt.Q_ARG("int", self._pl_freq_options.index(i)))
        self._pl_freq_callback(self.pl_freq)
        self._pl_freq_combo_box.currentIndexChanged.connect(
            lambda i: self.set_pl_freq(self._pl_freq_options[i]))
        # Create the radio buttons
        self.top_layout.addWidget(self._pl_freq_tool_bar)
        self._freq_tx_range = qtgui.Range(0, 915, .1, 145.99, 200)
        self._freq_tx_win = qtgui.RangeWidget(self._freq_tx_range, self.set_freq_tx, "'freq_tx'", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_layout.addWidget(self._freq_tx_win)
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
            ",".join(("", '')),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
            "",
        )
        self.uhd_usrp_sink_0.set_samp_rate(if_rate)
        self.uhd_usrp_sink_0.set_time_unknown_pps(uhd.time_spec(0))

        self.uhd_usrp_sink_0.set_center_freq(freq_tx*1e6, 0)
        self.uhd_usrp_sink_0.set_antenna("TX/RX", 0)
        self.uhd_usrp_sink_0.set_bandwidth(if_rate, 0)
        self.uhd_usrp_sink_0.set_gain(60, 0)
        self.qtgui_sink_x_0 = qtgui.sink_c(
            1024, #fftsize
            window.WIN_BLACKMAN_hARRIS, #wintype
            0, #fc
            if_rate, #bw
            "", #name
            True, #plotfreq
            True, #plotwaterfall
            True, #plottime
            True, #plotconst
            None # parent
        )
        self.qtgui_sink_x_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_win = sip.wrapinstance(self.qtgui_sink_x_0.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_0.enable_rf_freq(True)

        self.top_layout.addWidget(self._qtgui_sink_x_0_win)
        self.low_pass_filter_0_0 = filter.fir_filter_ccf(
            1,
            firdes.low_pass(
                1,
                if_rate,
                8000,
                2000,
                window.WIN_HAMMING,
                6.76))
        self.blocks_multiply_xx_0 = blocks.multiply_vcc(1)
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_ff(volume_gain)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_float*1, 'Z:\\schellberg\\indoor_exp_feb_2024\\SatelliteNetwork-main\\simpleTX_sim\\KJ5DKX_audio', False, 0, 0)
        self.blocks_file_source_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_add_xx_0 = blocks.add_vff(1)
        self.band_pass_filter_0 = filter.fir_filter_fff(
            1,
            firdes.band_pass(
                1,
                samp_rate,
                300,
                9000,
                200,
                window.WIN_HAMMING,
                6.76))
        self.analog_sig_source_x_0 = analog.sig_source_f(samp_rate, analog.GR_SIN_WAVE, pl_freq, .5, 0, 0)
        self.analog_nbfm_tx_0 = analog.nbfm_tx(
        	audio_rate=samp_rate,
        	quad_rate=if_rate,
        	tau=(75e-6),
        	max_dev=5e3,
        	fh=(-1),
                )
        self.analog_const_source_x_0 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, ptt)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_const_source_x_0, 0), (self.blocks_multiply_xx_0, 1))
        self.connect((self.analog_nbfm_tx_0, 0), (self.low_pass_filter_0_0, 0))
        self.connect((self.analog_sig_source_x_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.band_pass_filter_0, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.analog_nbfm_tx_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.band_pass_filter_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.blocks_multiply_xx_0, 0), (self.qtgui_sink_x_0, 0))
        self.connect((self.blocks_multiply_xx_0, 0), (self.uhd_usrp_sink_0, 0))
        self.connect((self.low_pass_filter_0_0, 0), (self.blocks_multiply_xx_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "NFM_xmt")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_usrp_rate(self):
        return self.usrp_rate

    def set_usrp_rate(self, usrp_rate):
        self.usrp_rate = usrp_rate
        self.set_if_rate(int(self.usrp_rate/3))

    def get_volume_gain(self):
        return self.volume_gain

    def set_volume_gain(self, volume_gain):
        self.volume_gain = volume_gain
        self.blocks_multiply_const_vxx_0.set_k(self.volume_gain)

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.analog_sig_source_x_0.set_sampling_freq(self.samp_rate)
        self.band_pass_filter_0.set_taps(firdes.band_pass(1, self.samp_rate, 300, 9000, 200, window.WIN_HAMMING, 6.76))

    def get_ptt(self):
        return self.ptt

    def set_ptt(self, ptt):
        self.ptt = ptt
        self.analog_const_source_x_0.set_offset(self.ptt)

    def get_pl_freq(self):
        return self.pl_freq

    def set_pl_freq(self, pl_freq):
        self.pl_freq = pl_freq
        self._pl_freq_callback(self.pl_freq)
        self.analog_sig_source_x_0.set_frequency(self.pl_freq)

    def get_if_rate(self):
        return self.if_rate

    def set_if_rate(self, if_rate):
        self.if_rate = if_rate
        self.low_pass_filter_0_0.set_taps(firdes.low_pass(1, self.if_rate, 8000, 2000, window.WIN_HAMMING, 6.76))
        self.qtgui_sink_x_0.set_frequency_range(0, self.if_rate)
        self.uhd_usrp_sink_0.set_samp_rate(self.if_rate)
        self.uhd_usrp_sink_0.set_bandwidth(self.if_rate, 0)

    def get_freq_tx(self):
        return self.freq_tx

    def set_freq_tx(self, freq_tx):
        self.freq_tx = freq_tx
        self.uhd_usrp_sink_0.set_center_freq(self.freq_tx*1e6, 0)




def main(top_block_cls=NFM_xmt, options=None):

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
