<html>
	<head>
		
	</head>
	<body>
	<br><br>
		<font face="arial" size="36">Done!
	<br><br>
		<a href="index.html"> New </a>
	<br><br>
	


<?php
$address = localhost;
$port = 4309;

$socket = socket_create(AF_INET, SOCK_STREAM, getprotobyname('tcp'));
socket_connect($socket, $address, $port);

$str = $_REQUEST['radius'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['thickness'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['angle'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['email'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['type'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['type1'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['area'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['Id'];
$message= $str."\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$str = $_REQUEST['lipzide'];
$message= "REMARK 547 ";
$str = $str."\n";
$message= $message.$str;
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);

$message= "exit\n";
$len = strlen($message);
$status = socket_sendto($socket, $message, $len, 0, $address, $port);
ini_set('max_execution_time', 300);
if($status !== FALSE)
{
    $message = '';
    $next = '';
    while ($next = socket_read($socket, 4096))
    {
        $message .= $next;
    }

    echo $message;
}
else
{
    echo "Failed";
}

socket_close($socket);
?>
</font>
</body>
<html>