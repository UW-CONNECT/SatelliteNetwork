#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Demo Spec An
# GNU Radio version: v3.8.2.0-57-gd71cd177

from distutils.version import StrictVersion

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print("Warning: failed to XInitThreads()")

from PyQt5 import Qt
from gnuradio import qtgui
from gnuradio.filter import firdes
import sip
from gnuradio import analog
from gnuradio import blocks
import pmt
from gnuradio import gr
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import zeromq
from gnuradio.qtgui import Range, RangeWidget
import math

from gnuradio import qtgui

class demo_spec_an(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Demo Spec An")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Demo Spec An")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
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

        self.settings = Qt.QSettings("GNU Radio", "demo_spec_an")

        try:
            if StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
                self.restoreGeometry(self.settings.value("geometry").toByteArray())
            else:
                self.restoreGeometry(self.settings.value("geometry"))
        except:
            pass

        ##################################################
        # Variables
        ##################################################
        self.snr_const_3d = snr_const_3d = -17
        self.snr = snr = 0
        self.samp_rate = samp_rate = 200000
        self.freq_offset = freq_offset = 0
        self.fft_size = fft_size = 1000

        ##################################################
        # Blocks
        ##################################################
        self._snr_range = Range(-100, 100, 1, 0, 200)
        self._snr_win = RangeWidget(self._snr_range, self.set_snr, 'Noise Level', "counter_slider", float)
        self.top_grid_layout.addWidget(self._snr_win)
        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, 'tcp://127.0.0.1:55555', 100, False, -1)
        self.qtgui_sink_x_0 = qtgui.sink_c(
            1024, #fftsize
            firdes.WIN_BLACKMAN_hARRIS, #wintype
            0, #fc
            samp_rate, #bw
            "", #name
            True, #plotfreq
            True, #plotwaterfall
            True, #plottime
            True #plotconst
        )
        self.qtgui_sink_x_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_win = sip.wrapinstance(self.qtgui_sink_x_0.pyqwidget(), Qt.QWidget)

        self.qtgui_sink_x_0.enable_rf_freq(False)

        self.top_grid_layout.addWidget(self._qtgui_sink_x_0_win)
        self.qtgui_number_sink_snr = qtgui.number_sink(
            gr.sizeof_float,
            0,
            qtgui.NUM_GRAPH_HORIZ,
            1
        )
        self.qtgui_number_sink_snr.set_update_time(0.10)
        self.qtgui_number_sink_snr.set_title("")

        labels = ['SNR', '', '', '', '',
            '', '', '', '', '']
        units = ['', '', '', '', '',
            '', '', '', '', '']
        colors = [("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"),
            ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black")]
        factor = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]

        for i in range(1):
            self.qtgui_number_sink_snr.set_min(i, -1)
            self.qtgui_number_sink_snr.set_max(i, 1)
            self.qtgui_number_sink_snr.set_color(i, colors[i][0], colors[i][1])
            if len(labels[i]) == 0:
                self.qtgui_number_sink_snr.set_label(i, "Data {0}".format(i))
            else:
                self.qtgui_number_sink_snr.set_label(i, labels[i])
            self.qtgui_number_sink_snr.set_unit(i, units[i])
            self.qtgui_number_sink_snr.set_factor(i, factor[i])

        self.qtgui_number_sink_snr.enable_autoscale(True)
        self._qtgui_number_sink_snr_win = sip.wrapinstance(self.qtgui_number_sink_snr.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_number_sink_snr_win)
        self.qtgui_number_sink_signal_amplitude = qtgui.number_sink(
            gr.sizeof_float,
            0,
            qtgui.NUM_GRAPH_HORIZ,
            1
        )
        self.qtgui_number_sink_signal_amplitude.set_update_time(0.10)
        self.qtgui_number_sink_signal_amplitude.set_title("")

        labels = ['Signal Amplitude', '', '', '', '',
            '', '', '', '', '']
        units = ['', '', '', '', '',
            '', '', '', '', '']
        colors = [("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"),
            ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black")]
        factor = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]

        for i in range(1):
            self.qtgui_number_sink_signal_amplitude.set_min(i, 0)
            self.qtgui_number_sink_signal_amplitude.set_max(i, 2)
            self.qtgui_number_sink_signal_amplitude.set_color(i, colors[i][0], colors[i][1])
            if len(labels[i]) == 0:
                self.qtgui_number_sink_signal_amplitude.set_label(i, "Data {0}".format(i))
            else:
                self.qtgui_number_sink_signal_amplitude.set_label(i, labels[i])
            self.qtgui_number_sink_signal_amplitude.set_unit(i, units[i])
            self.qtgui_number_sink_signal_amplitude.set_factor(i, factor[i])

        self.qtgui_number_sink_signal_amplitude.enable_autoscale(False)
        self._qtgui_number_sink_signal_amplitude_win = sip.wrapinstance(self.qtgui_number_sink_signal_amplitude.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_number_sink_signal_amplitude_win)
        self.qtgui_number_sink_0 = qtgui.number_sink(
            gr.sizeof_float,
            0,
            qtgui.NUM_GRAPH_HORIZ,
            1
        )
        self.qtgui_number_sink_0.set_update_time(0.10)
        self.qtgui_number_sink_0.set_title("")

        labels = ['Noise Amplitude', '', '', '', '',
            '', '', '', '', '']
        units = ['', '', '', '', '',
            '', '', '', '', '']
        colors = [("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"),
            ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black"), ("black", "black")]
        factor = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]

        for i in range(1):
            self.qtgui_number_sink_0.set_min(i, 0)
            self.qtgui_number_sink_0.set_max(i, 100)
            self.qtgui_number_sink_0.set_color(i, colors[i][0], colors[i][1])
            if len(labels[i]) == 0:
                self.qtgui_number_sink_0.set_label(i, "Data {0}".format(i))
            else:
                self.qtgui_number_sink_0.set_label(i, labels[i])
            self.qtgui_number_sink_0.set_unit(i, units[i])
            self.qtgui_number_sink_0.set_factor(i, factor[i])

        self.qtgui_number_sink_0.enable_autoscale(False)
        self._qtgui_number_sink_0_win = sip.wrapinstance(self.qtgui_number_sink_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_number_sink_0_win)
        self._freq_offset_range = Range(-5*samp_rate/fft_size, 5*samp_rate/fft_size, 1, 0, 200)
        self._freq_offset_win = RangeWidget(self._freq_offset_range, self.set_freq_offset, 'Frequency Offset', "counter_slider", float)
        self.top_grid_layout.addWidget(self._freq_offset_win)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_rms_xx_0_0 = blocks.rms_cf(.0001)
        self.blocks_rms_xx_0 = blocks.rms_cf(.0001)
        self.blocks_nlog10_ff_0 = blocks.nlog10_ff(20, 1, 0)
        self.blocks_multiply_xx_0 = blocks.multiply_vcc(1)
        self.blocks_file_source_1 = blocks.file_source(gr.sizeof_gr_complex*1, 'C:\\Users\\schellberg\\Documents\\schellberg\\SatelliteNetwork\\simpleTX_sim\\thousand5_tests\\SF9_noOnes', True, 0, 0)
        self.blocks_file_source_1.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_source_0 = blocks.file_source(gr.sizeof_gr_complex*1, 'C:\\Users\\schellberg\\Documents\\schellberg\\SatelliteNetwork\\simpleTX_sim\\doppler_sim_scripts\\GNURADIO_DOP_REF', True, 0, 0)
        self.blocks_file_source_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_divide_xx_0 = blocks.divide_ff(1)
        self.blocks_add_xx_0 = blocks.add_vcc(1)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, math.pow(10.0,snr/20.0), 42)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_rms_xx_0, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.blocks_divide_xx_0, 0), (self.blocks_nlog10_ff_0, 0))
        self.connect((self.blocks_file_source_0, 0), (self.blocks_multiply_xx_0, 1))
        self.connect((self.blocks_file_source_1, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.blocks_file_source_1, 0), (self.blocks_rms_xx_0_0, 0))
        self.connect((self.blocks_multiply_xx_0, 0), (self.qtgui_sink_x_0, 0))
        self.connect((self.blocks_multiply_xx_0, 0), (self.zeromq_pub_sink_0, 0))
        self.connect((self.blocks_nlog10_ff_0, 0), (self.qtgui_number_sink_snr, 0))
        self.connect((self.blocks_rms_xx_0, 0), (self.blocks_divide_xx_0, 1))
        self.connect((self.blocks_rms_xx_0, 0), (self.qtgui_number_sink_0, 0))
        self.connect((self.blocks_rms_xx_0_0, 0), (self.blocks_divide_xx_0, 0))
        self.connect((self.blocks_rms_xx_0_0, 0), (self.qtgui_number_sink_signal_amplitude, 0))
        self.connect((self.blocks_throttle_0, 0), (self.blocks_multiply_xx_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "demo_spec_an")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_snr_const_3d(self):
        return self.snr_const_3d

    def set_snr_const_3d(self, snr_const_3d):
        self.snr_const_3d = snr_const_3d

    def get_snr(self):
        return self.snr

    def set_snr(self, snr):
        self.snr = snr
        self.analog_noise_source_x_0.set_amplitude(math.pow(10.0,self.snr/20.0))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)
        self.blocks_throttle_0_0.set_sample_rate(self.samp_rate)
        self.qtgui_sink_x_0.set_frequency_range(0, self.samp_rate)

    def get_freq_offset(self):
        return self.freq_offset

    def set_freq_offset(self, freq_offset):
        self.freq_offset = freq_offset

    def get_fft_size(self):
        return self.fft_size

    def set_fft_size(self, fft_size):
        self.fft_size = fft_size





def main(top_block_cls=demo_spec_an, options=None):

    if StrictVersion("4.5.0") <= StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    def quitting():
        tb.stop()
        tb.wait()

    qapp.aboutToQuit.connect(quitting)
    qapp.exec_()

if __name__ == '__main__':
    main()
