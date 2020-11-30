function message_REX(message)
global udpCom
pnet(udpCom.sock, 'write', [message '>> >>']);
pnet(udpCom.sock, 'writepacket', udpCom.rexip, udpCom.port);
