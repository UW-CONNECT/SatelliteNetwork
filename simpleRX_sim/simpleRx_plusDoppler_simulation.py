#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# GNU Radio version: 3.10.8.0

from PyQt5 import Qt
from gnuradio import qtgui
from gnuradio import analog
from gnuradio import blocks
import pmt
from gnuradio import channels
from gnuradio.filter import firdes
from gnuradio import gr
from gnuradio.fft import window
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import zeromq
import sip



class simpleRx_plusDoppler_simulation(gr.top_block, Qt.QWidget):

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

        self.settings = Qt.QSettings("GNU Radio", "simpleRx_plusDoppler_simulation")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 200e3
        self.freq_slope = freq_slope = 50
        self.sweep_time_seconds = sweep_time_seconds = (samp_rate/2)/freq_slope

        ##################################################
        # Blocks
        ##################################################

        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, 'tcp://127.0.0.1:55555', 100, False, (-1), '', True, True)
        self.qtgui_sink_x_1 = qtgui.sink_c(
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
        self.qtgui_sink_x_1.set_update_time(1.0/10)
        self._qtgui_sink_x_1_win = sip.wrapinstance(self.qtgui_sink_x_1.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_1.enable_rf_freq(False)

        self.top_layout.addWidget(self._qtgui_sink_x_1_win)
        self.channels_channel_model_1 = channels.channel_model(
            noise_voltage=.003,
            frequency_offset=0,
            epsilon=1.0,
            taps=[1],
            noise_seed=0,
            block_tags=False)
        self.blocks_vco_c_0 = blocks.vco_c(samp_rate, (((2*3.14)*freq_slope*sweep_time_seconds) ), 1)
        self.blocks_throttle2_0 = blocks.throttle( gr.sizeof_float*1, samp_rate, True, 0 if "auto" == "auto" else max( int(float(0.1) * samp_rate) if "auto" == "time" else int(0.1), 1) )
        self.blocks_multiply_xx_0 = blocks.multiply_vcc(1)
        self.blocks_multiply_const_vxx_1 = blocks.multiply_const_cc(.005)
        self.blocks_file_source_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, 'J:\\schellberg\\indoor_exp_feb_2024\\TEST_ISOLATED\\0HzS_SF_7N_128BW_2500FS_200000NPKTS_50PLEN_100CR_0\\trial1', False, 0, 0)
        self.blocks_file_source_0_0.set_begin_tag(pmt.PMT_NIL)
        self.analog_sig_source_x_0 = analog.sig_source_f(samp_rate, analog.GR_SAW_WAVE, (1/sweep_time_seconds), 1, (-1/2), 0)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_sig_source_x_0, 0), (self.blocks_throttle2_0, 0))
        self.connect((self.blocks_file_source_0_0, 0), (self.blocks_multiply_const_vxx_1, 0))
        self.connect((self.blocks_multiply_const_vxx_1, 0), (self.blocks_multiply_xx_0, 0))
        self.connect((self.blocks_multiply_xx_0, 0), (self.channels_channel_model_1, 0))
        self.connect((self.blocks_throttle2_0, 0), (self.blocks_vco_c_0, 0))
        self.connect((self.blocks_vco_c_0, 0), (self.blocks_multiply_xx_0, 1))
        self.connect((self.channels_channel_model_1, 0), (self.qtgui_sink_x_1, 0))
        self.connect((self.channels_channel_model_1, 0), (self.zeromq_pub_sink_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "simpleRx_plusDoppler_simulation")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_sweep_time_seconds((self.samp_rate/2)/self.freq_slope)
        self.analog_sig_source_x_0.set_sampling_freq(self.samp_rate)
        self.blocks_throttle2_0.set_sample_rate(self.samp_rate)
        self.qtgui_sink_x_1.set_frequency_range(0, self.samp_rate)

    def get_freq_slope(self):
        return self.freq_slope

    def set_freq_slope(self, freq_slope):
        self.freq_slope = freq_slope
        self.set_sweep_time_seconds((self.samp_rate/2)/self.freq_slope)

    def get_sweep_time_seconds(self):
        return self.sweep_time_seconds

    def set_sweep_time_seconds(self, sweep_time_seconds):
        self.sweep_time_seconds = sweep_time_seconds
        self.analog_sig_source_x_0.set_frequency((1/self.sweep_time_seconds))




def main(top_block_cls=simpleRx_plusDoppler_simulation, options=None):

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
