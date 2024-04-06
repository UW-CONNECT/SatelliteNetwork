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
        self.samp_rate = samp_rate = 200000

        ##################################################
        # Blocks
        ##################################################

        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 1, 'tcp://127.0.0.1:55555', 100, False, (-1), '', True, True)
        self.qtgui_sink_x_0 = qtgui.sink_c(
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
        self.qtgui_sink_x_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_win = sip.wrapinstance(self.qtgui_sink_x_0.qwidget(), Qt.QWidget)

        self.qtgui_sink_x_0.enable_rf_freq(False)

        self.top_layout.addWidget(self._qtgui_sink_x_0_win)
        self.channels_channel_model_0 = channels.channel_model(
            noise_voltage=0.025,
            frequency_offset=0.0,
            epsilon=1.0,
            taps=[1.0 + 1.0j],
            noise_seed=0,
            block_tags=False)
        self.blocks_throttle2_0 = blocks.throttle( gr.sizeof_gr_complex*1, samp_rate, True, 0 if "auto" == "auto" else max( int(float(0.1) * samp_rate) if "auto" == "time" else int(0.1), 1) )
        self.blocks_multiply_xx_2 = blocks.multiply_vcc(1)
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_cc(.008)
        self.blocks_file_source_3 = blocks.file_source(gr.sizeof_gr_complex*1, 'J:\\schellberg\\indoor_exp_feb_2024\\doppler_simulation_files\\GNURADIO_linear_8k_sep', True, 0, 0)
        self.blocks_file_source_3.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_source_1 = blocks.file_source(gr.sizeof_gr_complex*1, 'J:\\schellberg\\indoor_exp_feb_2024\\experiment_data\\SF_7N_128BW_40000FS_200000NPKTS_100PLEN_100CR_0\\trial1', False, 0, 0)
        self.blocks_file_source_1.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_sink_1 = blocks.file_sink(gr.sizeof_gr_complex*1, 'C:\\Users\\schellberg\\Documents\\schellberg\\Standard_LoRa\\sf7_with_doppler', False)
        self.blocks_file_sink_1.set_unbuffered(False)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_1, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.blocks_file_source_3, 0), (self.blocks_multiply_xx_2, 1))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.blocks_multiply_xx_2, 0))
        self.connect((self.blocks_multiply_xx_2, 0), (self.blocks_throttle2_0, 0))
        self.connect((self.blocks_throttle2_0, 0), (self.channels_channel_model_0, 0))
        self.connect((self.channels_channel_model_0, 0), (self.blocks_file_sink_1, 0))
        self.connect((self.channels_channel_model_0, 0), (self.qtgui_sink_x_0, 0))
        self.connect((self.channels_channel_model_0, 0), (self.zeromq_pub_sink_0, 0))


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
        self.blocks_throttle2_0.set_sample_rate(self.samp_rate)
        self.qtgui_sink_x_0.set_frequency_range(0, self.samp_rate)




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
